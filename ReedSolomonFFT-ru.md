[Помехо-защищённое кодирование][ECC-ru] (ECC - [error-correcting codes][ECC]) используется в протоколах передачи данных (от древних аналоговых модемов до 4G LTE и цифрового ТВ), хранении данных ([DVD], HDD, SSD, [RAID-6]), наконец в прикладных программах, реализующих в том или ином виде "recovery record" ([RAR], [PAR2], [RSC32]).

Наиболее известными в этой области являются [коды Рида-Соломона][Reed-Solomon codes-ru] ([Reed-Solomon codes]). Я расскажу о том, как они устроены, зачем тут нужны поля Галуа и как использовать FFT для RS-кодирования за время O(n*log(n)), что позволяет добиться скорости в сотни мегабайт/сек даже при кодировании групп, состоящих из многих миллионов блоков. Вам потребуются лишь минимальные познания в математике: [полиномиальная арифметика][polynomial arithmetic], [умножение][matrix multiplication] и [инвертирование][matrix inverse] матриц, [матрица Вандермонда][Vandermonde matrix].

* [Общие принципы ECC](#general)
* [Коды Рида-Соломона](#reed-solomon)
  * [Полиномиальный подход](#polynomial)
  * [Матричный подход](#matrix)
  * [Матричный подход на базе матриц Вандермонда](#vandermonde)
* [Быстрые алгоритмы](#fast)
  * [DFT](#dft)
  * [FFT и быстрая полиномиальная арифметика](#fft)
  * [Кодирование/декодирование за время O(n*log<sup>2</sup>(n))](#simple)
  * [Кодирование за время O(n*log(n))](#encode)
  * [Декодирование за время O(n*log(n))](#decode)
* [Поля Галуа](#galua)


<a name="general"/>

# Общие принципы ECC

Идея любого помехо-защищённого кода проста - имея k подлежащих защите элементов (слов) данных, мы вычисляем по ним некие n слов (n>k), передаём их (или храним в течении некоторого времени), и при восстановлении вычисляем k исходных слов назад, даже если некоторые из этих n переданных/сохранённых слов были утеряны или искажены. Это обозначается как (n,k) код. Из теории информации очевидно, что для гарантированного восстановления необходимо получить неискажёнными хотя бы k слов, и что имеют смысл только коды с n>k.

Теперь немного терминологии.

Искажения при передаче данных делятся на два типа:
1) ошибки (errors): это неправильно переданные данные, про которые получателю неизвестно, корректны они или нет
2) стирания (erasures): это слова, которые точно были переданы неправильно, и потому игнорируются в процессе декодирования

Традиционно ECC были более нацелены на восстановление ошибок (errors), однако по мере роста размеров групп кодируемых данных становилось всё более выгодно разбивать их на блоки, сохранять контрольную сумму каждого блока и при её несовпадении объявлять все слова в блоке потерянными (т.е. erasures). Причина тут в том, что каждое дополнительное слово ECC позволяет восстановить одно стирание (т.е. неправильное слово, когда мы точно знаем что оно неправильно), но для исправления каждой ошибки необходимо два слова ECC. В целом, когда полученные данные содержат E ошибок и W стираний, для восстановления необходимо 2E+W слов ECC.


Все коды делятся на две группы:
1) максимально-разделимые (MDS, maximum distance separable): это коды, гарантирующие восстановление данных при n-k>=2E+W. В частности, имея дело только с erasures, они гарантируют восстановление при n-k>=W, т.е. если удалось передать без искажений любые k из n слов кода (и мы знаем, какие слова были переданы верно). Понятно, что это наиболее "качественные" коды, позволяющие восстановить теоретически возможный максимум.
2) non-MDS коды: все прочие. В последнее время многие семейства таких кодов (LDPC, Raptor, fountain...) приобрели бешенную популярность из-за своих характеристик - во-первых, число лишних блоков, требуемых для восстановления данных, может достигать всего 1-2 на миллион блоков (что делает на практике разницу с MDS-кодами несущественной), во-вторых скорость кодирования на порядки превышает скорость классических O(n<sup>2</sup>) алгоритмов для кодов Рида-Соломона, в-третьих у из-за различных дополнительных возможностей (soft-декодирования, фонтанного кодирования...). Как видите, данная работа как раз описывает как создать коды Рида-Соломона имеющие скорость O(n*log(n)), что почти нивелирует преимущество non-MDS кодов в скорости, но при этом всё ещё являющиеся стопроцентными MDS.


С другой стороны, коды можно разделить на
1) систематические: те, где в состав n передаваемых слов входят k исходных слов. В таких кодах остальные кодовые слова называются словами чётности (parity words). Некоторые приложения (скажем добавление ECC к существующим файлам на диске) требуют применения именно таких кодов. Если исходные данные не повреждены, такие коды позволяют обойтись вообще без декодирования.
2) несистематические: передаваемые данные не включают в себя исходные. Декодирование всегда необходимо. Само по себе, очевидно, это не даёт никаких преимуществ перед систематическими кодами. Их используют только в тех случаях, когда они оказываются вычислительно дешевле аналогичного систематического кода, и при этом приложение не требует, чтобы код был систематическим.


<a name="reed-solomon"/>

# Коды Рида-Соломона

<a name="polynomial"/>

## Полиномиальный подход

Пререквизиты: [сложение/умножение/деление полиномов][polynomial arithmetic]

Наиболее известными из MDS-кодов являются коды Рида-Соломона (Reed-Solomon codes), являющиеся частным случаем BCH-кодов. Существует множество алгоритмов, относящихся к этой группе, но видимо их все можно разделить на два класса:
1) сводящиеся к умножению и делению полиномов
2) сводящиеся к умножению и "делению" матриц

Первая группа, как я понимаю, является классическими алгоритмами, поэтому мы для начала кратко рассмотрим её. Вы можете пропустить этот раздел, если вас интересует только алгоритм FastECC.

Назовём полиномами "длины" n все полиномы, которые можно представить в виде p(x)=a<sub>0</sub>+a<sub>1</sub>x+...+a<sub>n-1</sub>x<sup>n-1</sup>, где любые коэффициенты a<sub>i</sub> могут быть нулевыми. Очевидно, что это множество включает все полиномы степеней от 0 до n-1, и существует взаимно-однозначное соответствие между множеством полиномов "длины" n и множеством векторов длины n.

Первый алгоритм "полиномиального" RS-кодирования очень прост: рассмотрим k входных слов как полином `a` "длины" k и умножим его на специальный полином `g` степени `n-k`, который мы назовём генератором. Произведение `a*g` является полиномом "длины" n - это и есть передаваемое сообщение. При этом разработанные алгоритмы RS-декодирования позволяют восстановить исходный полином `a` при условии n-k>=2E+W, что делает его MDS-кодом. Но, как видите, передаваемые данные (коэффициенты полинома `a*g`) не содержат исходные данные (коэффициенты полинома `a`), т.е. этот код - несистематический.

Аналогичный систематический код выглядит следующим образом: умножим `a` на x<sup>n-k</sup>, найдём остаток от деления произведения на `g` и вычтем его из произведения: <code>p = a*x<sup>n-k</sup>, r = p rem g, q = p-r</code>. Полином `q` одновременно обладает двумя замечательными свойствами - он делится на g и при этом k его старших коэффициентов совпадают с коэффициентами полинома `a`, т.е. входными данными. Второе свойство делает полином q систематическим кодом. Первое же свойство означает, что к нему применимы точно те же алгоритмы декодирования, что и к предыдущему (несистематическому) коду - если полученный этими алгоритмами результат умножить на g и взять k старших коэффициентов произведения, то мы получим результат декодирования систематического кода!

Как видите при кодировании у нас есть три time-consuming задачи: генерация полинома g, умножение полиномов и деление полиномов. Умножение полиномов можно выполнить за время `O(n*log(n))`, где n - степень результата. Генерация полинома g требует <code>O(n&#42;log<sup>2</sup>(n))</code> времени, однако выполняется только однажды и при использовании блоков, содержащих минимум сотни элементов (т.е. сотни байт данных), оказывается куда меньше времени умножения/деления. Наконец, деление тоже требует <code>O(n&#42;log<sup>2</sup>(n))</code> времени, однако ЕМНИП, основное время при этом тратится на генерацию "обратного" полинома, после чего происходит умножение на него за время `O(n*log(n))`. Поскольку нам нужно только деление на фиксированный полином g, "обратный" к нему полином нужно сгенерить лишь однократно, и время этой операции опять же оказывается существенно меньше времени, требуемого на дальнейшее умножение сотен полиномов, соответствующих элементам блоков данных.

В результате, с практической точки зрения, даже алгоритм систематического кодирования работает за время `O(n*log(n))`. Программа [RSC32], которая судя по всему реализует именно его, достигает скорости в сотни мегабайт в секунду, и порядка 100 МБ/с при декодировании.


<a name="matrix"/>

## Матричный подход

Пререквизиты: [умножение матрицы на вектор/матрицу][matrix multiplication], [обратная матрица][matrix inverse]

Матричный подход к кодам Рида-Соломона сводится к следующему алгоритму: возьмём вектор A длины k, представляющий входные данные, и умножим его на заранее заданную матрицу G размером `k*n`. Полученный вектор длины n и будет представлять кодовое слово. Легко видеть, что это систематический код, если G включает в себя единичную подматрицу размера `k*k`, и такой код можно описать оставшейся частью матрицы размером `k*(n-k)` для вычисления слов чётности.

Алгоритм декодирования (при наличии только erasures) прост и элегантен: возьмём вектор B - любые k из успешно полученных значений, и выпишем E - подматрицу G размером `k*k`, на которую нужно умножить A, чтобы их получить: `B=E*A`. Для восстановления A достаточно умножить B на обратную к E матрицу: <code>A=B*E<sup>-1</sup></code>.

Нетрудно видеть, что однозначное декодирование по k сохранившимся элементам возможно тогда и только тогда, когда генерирующая их подматрица G обратима, а MDS-свойством обладают те и только те коды, у которых любые k строк матрицы G составляют обратимую матрицу.

AFAIK, любая [матрица Коши][Cauchy matrix] размером `k*(n-k)` даёт систематический MDS-код. AFAIK, этот код используется в RAR5.


<a name="vandermonde"/>

## Матричный подход на базе матриц Вандермонда

Пререквизиты: [матрица Вандермонда][Vandermonde matrix]

Значение полинома "длины" k (см. выше) в точке x можно вычислить по формуле p(x)=a<sub>0</sub>+a<sub>1</sub>*x+...a<sub>k-1</sub>*x<sup>k-1</sup>, что эквивалентно скалярному произведению векторов A=(a<sub>0</sub>,a<sub>1</sub>...a<sub>k-1</sub>) и (1,x...x<sup>k-1</sup>). Вычислить значения полинома сразу в n точках (x<sub>0</sub>,x<sub>1</sub>...x<sub>n-1</sub>) можно, умножив построенную из них матрицу Вандермонда на тот же вектор A:

<pre>
P = V*A, где
P = (p(x<sub>0</sub>), p(x<sub>1</sub>) ... p(x<sub>n-1</sub>))
V = Vandermonde<sub>k</sub>(x<sub>0</sub>, x<sub>1</sub> ... x<sub>n-1</sub>)
A = (a<sub>0</sub>, a<sub>1</sub> ... a<sub>k-1</sub>)
</pre>

Нетрудно сообразить, что в обратную сторону - вычислить коэффициенты полинома "длины" k по его значениям в k различных точках, можно, умножив P на матрицу, обратную к V. Матрица Вандермонда всегда обратима (если x<sub>i</sub> попарно различны), следовательно коэффициенты полинома "длины" k всегда и однозначно восстановимы по его значениям в любых k точках.

Отсюда следует что существует взаимно-однозначное соответствие между двумя равномощными множествами: множеством векторов длины k, рассматриваемыми как коэффициенты полиномов "длины" k, и множеством векторов той же длины k, рассматриваемыми как значения полиномов "длины" k в некоторых k фиксированных различных точках.

Любая подматрица `k*k` матрицы Вандермонда `k*n` сама по себе является матрицей Вандермонда и следовательно обратима.

Теперь мы можем рассмотреть новые классы кодов Рида-Соломона, описываемые соответствующими генерирующими матрицами G:
- матрица Вандермонда размером `k*n` даёт несистематический MDS-код. Такой код эквивалентен тому, что исходные данные трактуются как коэффициенты полинома "длины" k, а кодовое слово представляет собой значения этого полинома в некоторых n точках. Декодирование эквивалентно восстановлению коэффициентов полинома "длины" k по его значениям в k точках. Это всегда возможно, следовательно мы имеем дело с MDS-кодом.
- матрица Вандермонда размером `k*(n-k)` даёт систематический код. Он эквивалентен включению в состав кодового слова k коэффициентов полинома и его значений в n-k точках, и попытке восстановления по любым k из этих n значений. Этот код был рекомендован Планком в его знаменитом тюториале 1996-го года и реализован в формате [PAR2]. К сожалению, позже выяснилось, что он не всегда является MDS-кодом, хотя видимо достаточно близок к этому с практической точки зрения.
- и наконец систематический MDS-код, реализованный в FastECC: используем исходные k слов как значения некоего полинома "длины" k в k точках, и вычислим слова чётности как его значения в других n-k точках. В результате, кодовое слово будет состоять из значений полинома длины k в n точках, и очевидно, что восстановление возможно по любым k сохранившимся значениям. При этом для кодирования нужно сначала умножить исходные данные на обратную матрицу Вандермонда размера `k*k` (чтобы получить коэф. полинома), и затем умножить результат на матрицу Вандермонда размера `k*(n-k)`. На практике, произведение этих двух матриц лучше вычислить один раз заранее, что и даёт нам матричный код с генерирующей матрицей <code>G = Vandermonde<sub>k</sub><sup>-1</sup>(x<sub>0</sub>,x<sub>1</sub>...x<sub>k-1</sub>) * Vandermonde<sub>k</sub>(x<sub>k</sub>,x<sub>k+1</sub>...x<sub>n-1</sub>)</code>. При декодировании мы вычисляем коэф. полинома по сохранившимся значениям, а по ним вычисляем значения в k точках, соответствующих исходным данным.

Осталось рассмотреть вопрос о скорости кодирования/декодирования. Кодирование исходных данных размером k элементов сводится к умножению на матрицу G размером `k*n`, что требует `O(k*n)` арифметических операций. При фиксированном соотношении k/n получается, что кодирование n элементов требует O(n<sup>2</sup>) времени, т.е. O(n) на один элемент. А это означает, что скорость кодирования обратно пропорциональна значениям n и k - скажем, код с (n,k)=(1000,500) будет вдесятеро медленней кода (100,50). На нынешних процессорах кодирование с n=1000 работает со скоростью в десятки МБ/с, т.е. уже медленнее, чем сами HDD. Кодирование с десятками тысяч и тем более миллионами блоков становится ещё более непрактичным. И вот тут-то нам на выручку и приходит FFT, к рассмотрению которого мы наконец и переходим.


<a name="fast"/>

# Быстрые алгоритмы

<a name="dft"/>

## DFT

Начнём с рассмотрения DFT (ДПФ - дискретное преобразование Фурье). Это математическая формула, которой можно найти десятки применений и соответственно толкований, нас же интересует следующее: DFT<sub>N</sub>(A,r) - это значения полинома, представленного вектором коэффициентов A длины N, в N точках (r<sup>0</sup>,r<sup>1</sup>...r<sup>n-1</sup>). При этом r должен быть [первообразным корнем N-й степени из 1][primitive root of unity] (такие корни есть, например, в поле комплексных чисел и в некоторых полях Галуа), что гарантирует отсутствие совпадений среди этих точек.

Таким образом, DFT N-го порядка преобразует вектор длины N в другой вектор той же длины, и эквивалентно умножению на матрицу Вандермонда специального вида:

<pre>
DFT<sub>N</sub>(A,r) = V*A, где
V = Vandermonde<sub>N</sub>(r<sup>0</sup>, r<sup>1</sup> ... r<sup>N-1</sup>)
</pre>

Поскольку все точки x<sub>i</sub> = r<sup>i</sup> различны, DFT является взаимно однозначным соответствием и всегда обратимо.

Обратное ДПФ (Inverse DFT), а именно IDFT<sub>N</sub>(A,r) собственно и выполняет обратное преобразование - из значений полинома в точках r<sup>i</sup> в его коэффициенты. Самое удивительное, что оба преобразования одинаковы с точностью до скалярного коэффициента: <code>IDFT<sub>N</sub>(A,r) = DFT<sub>N</sub>(A,r<sup>-1</sup>) / N</code>! Поэтому на базе любой функции, вычисляющей DFT, легко конструируется функция, вычисляющая IDFT.


<a name="fft"/>

## FFT и быстрая полиномиальная арифметика

Если DFT - это формула, то FFT (БПФ - быстрое преобразование Фурье) - это целый класс алгоритмов, обеспечивающих её вычисление за время `O(N*log(N))` в различных частных случаях. Аналогично, Inverse FFT (обратное БПФ) - это алгоритм, вычисляющий IDFT за время `O(N*log(N))`, и такие алгоритмы тривиально конструируются из алгоритмов FFT добавлением скалярного деления результата на N.

На базе FFT строятся быстрые алгоритмы полиномиальной арифметики:
* [умножение полиномов][fast polynomial multiplication] за время `O(n*log(n))`
* [деление полиномов][fast polynomial division] за время <code>O(n&#42;log<sup>2</sup>(n))</code>
* [polynomial evaluation][fast polynomial evaluation] за время <code>O(n&#42;log<sup>2</sup>(n))</code>: вычисление значений полинома "длины" n в n точках
* [polynomial interpolation][fast polynomial interpolation] за время <code>O(n&#42;log<sup>2</sup>(n))</code>: вычисление коэффициентов полинома "длины" n по его значениям в n точках

По приведённым ссылкам вы найдёте различные университетские лекции, описывающие эти алгоритмы. Обычным же образом, без использования FFT, все эти операции выполняются за время O(n<sup>2</sup>).

Дальнейшее чтение:
- [онлайн-описание][fft-online] основных алгоритмов FFT, за исключением [PFA]


<a name="simple"/>

## Кодирование/декодирование за время O(n*log<sup>2</sup>(n))

Имея в своём распоряжении эти быстрые polynomial evaluation/interpolation алгоритмы, легко реализовать кодирование и декодирование за время <code>O(n&#42;log<sup>2</sup>(n))</code>.

**Кодирование**: будем считать k входных слов значениями некоего полинома "длины" k в k определённых точках. Найдём коэффициенты этого полинома с помощью fast polynomial interpolation, и затем найдём его значения в других n-k точках с помощью fast polynomial evaluation. Обе операции укладыватся в <code>O(n&#42;log<sup>2</sup>(n))</code>, ЧТД.

**Декодирование**: почти не отличается от кодирования. Мы имеем k значений полинома в неких точках и знаем, что его степень меньше k. Это позволяет нам вычислить коэффициенты этого полинома с помощью fast polynomial interpolation, и затем найти значения в тех точках, которые соответствуют исходным данным, и тем самым восстановить их.

Как видите, эти алгоритмы легко описать и очень просто реализовать, если в вашем распоряжении есть библиотека быстрой полиномиальной арифметики в полях Галуа (а такие библиотеки можно найти даже для GPU). При этом можно рассчитывать на скорость порядка 100 МБ/с на CPU и 1 ГБ/c на GPU (при кодировании миллиона блоков).


<a name="encode"/>

## Кодирование за время O(n*log(n))

Таким образом, канва предлагаемого алгоритма кодирования становится понятна - трактуя входные данные A как k значений полинома "длины" k, вычислим B = IDFT<sub>k</sub>(A,r) и получим его коэффициенты. Дополнив эти коэффициенты старшими нулями, расширим B до размера n, и вычислим DFT<sub>n</sub>(B,r). Входные k значений и n-k из выходных значений являются значениями полинома степени k в n различных точках, и следовательно позволяют восстановить его по любым k сохранившимся значениям. Таким образом, с помощью FFT/IFFT мы можем осуществить систематическое MDS-кодирование по схеме FastECC за время `O(n*log(n))`. Алгоритм декодирования за время `O(n*log(n))` несколько сложнее, но сначала мы более детально рассмотрим алгоритм кодирования.


<a name="decode"/>

## Декодирование за время O(n*log(n))

Описание пока не готово. На английском описано в [README#how](README.md#how) и последующих секциях.


<a name="galua"/>

# Поля Галуа

Подробней: [поля Галуа][Galua field]

Как понимаете, всё это хорошо звучит абстрактно-математически, на практике же получится следующее: пусть нам надо закодировать N байт, так что будем их считать 8-битными числами. Тогда дополнительные M значений, вычисленные по любой из этих схем, могут выйти за пределы диапазона 0..255, и нам придётся хранить скажем `2*M` дополнительных байт, и соответственно для успешного восстановления данных тоже может потребоваться больше чем N байт. В последней же схеме дополнительные значения могут вообще получиться иррациональными. Даже если рассматривать входные данные как значения float/double, всё равно для однозначного восстановления по N значениям доп. данные потребуется хранить с большей точностью, чем входные.

И вот тут-то нам на помощь и приходят поля Галуа (GF - Galois Field). Детали вы найдёте в Википедии, сейчас достаточно сказать что поле Галуа GF(p<sup>n</sup>), где p - простое число, а n>=1 - любое натуральное, состоит из чисел от 0 до p<sup>n</sup>-1 и определяет операции сложения, вычитания, умножения и деления так, что они определены для любой пары элементов этого поля (за исключением деления на 0), однозначны и дают в результате элемент того же поля (т.е. значение от 0 до p<sup>n</sup>)-1.

А это уже означает, что если в любой нашей схеме кодирования N входных значений принадлежат некоторому полю GF(p<sup>n</sup>), то и выходные значения принадлежат тому же полю - т.е. `k*N` входных бит закодированы в `k*(N+M)` выходных бит, и любых сохранившихся N k-битных значений из них достаточно для восстановления входных данных. На практике обычно используют поля вида GF(2<sup>n</sup>), поскольку их элементы - обычные n-битные слова. Однако алгоритмы FFT с такими полями не дружат, поэтому для O(n*log(n)) кодирования мы будем использовать поля GF(p) - т.е. с n=1, которые требуют несложного перекодирования для представления произвольных бинарных данных.


[ECC-ru]: https://ru.wikipedia.org/wiki/%D0%9E%D0%B1%D0%BD%D0%B0%D1%80%D1%83%D0%B6%D0%B5%D0%BD%D0%B8%D0%B5_%D0%B8_%D0%B8%D1%81%D0%BF%D1%80%D0%B0%D0%B2%D0%BB%D0%B5%D0%BD%D0%B8%D0%B5_%D0%BE%D1%88%D0%B8%D0%B1%D0%BE%D0%BA#.D0.9A.D0.BE.D0.B4.D1.8B_.D0.BE.D0.B1.D0.BD.D0.B0.D1.80.D1.83.D0.B6.D0.B5.D0.BD.D0.B8.D1.8F_.D0.B8_.D0.B8.D1.81.D0.BF.D1.80.D0.B0.D0.B2.D0.BB.D0.B5.D0.BD.D0.B8.D1.8F_.D0.BE.D1.88.D0.B8.D0.B1.D0.BE.D0.BA
[ECC]: https://en.wikipedia.org/wiki/Error_detection_and_correction#Error-correcting_codes
[Reed-Solomon codes-ru]: https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D0%B4_%D0%A0%D0%B8%D0%B4%D0%B0_%E2%80%94_%D0%A1%D0%BE%D0%BB%D0%BE%D0%BC%D0%BE%D0%BD%D0%B0
[Reed-Solomon codes]: https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction

[DVD]: https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D0%B4_%D0%A0%D0%B8%D0%B4%D0%B0_%E2%80%94_%D0%A1%D0%BE%D0%BB%D0%BE%D0%BC%D0%BE%D0%BD%D0%B0#.D0.97.D0.B0.D0.BF.D0.B8.D1.81.D1.8C_.D0.BD.D0.B0_CD-ROM
[RAID-6]: https://en.wikipedia.org/wiki/Standard_RAID_levels#RAID_6
[RAR]: https://en.wikipedia.org/wiki/RAR_(file_format)
[PAR2]: https://en.wikipedia.org/wiki/Parchive
[RSC32]: https://www.livebusinesschat.com/smf/index.php?board=399.0

[polynomial arithmetic]: https://en.wikipedia.org/wiki/Polynomial#Arithmetic
[matrix multiplication]: https://ru.wikipedia.org/wiki/%D0%A3%D0%BC%D0%BD%D0%BE%D0%B6%D0%B5%D0%BD%D0%B8%D0%B5_%D0%BC%D0%B0%D1%82%D1%80%D0%B8%D1%86
[matrix inverse]: https://ru.wikipedia.org/wiki/%D0%9E%D0%B1%D1%80%D0%B0%D1%82%D0%BD%D0%B0%D1%8F_%D0%BC%D0%B0%D1%82%D1%80%D0%B8%D1%86%D0%B0

[Cauchy matrix]: https://ru.wikipedia.org/wiki/%D0%9C%D0%B0%D1%82%D1%80%D0%B8%D1%86%D0%B0_%D0%9A%D0%BE%D1%88%D0%B8_(%D0%BB%D0%B8%D0%BD%D0%B5%D0%B9%D0%BD%D0%B0%D1%8F_%D0%B0%D0%BB%D0%B3%D0%B5%D0%B1%D1%80%D0%B0)
[Vandermonde matrix]: https://en.wikipedia.org/wiki/Vandermonde_matrix
[primitive root of unity]: https://ru.wikipedia.org/wiki/%D0%9F%D0%B5%D1%80%D0%B2%D0%BE%D0%BE%D0%B1%D1%80%D0%B0%D0%B7%D0%BD%D1%8B%D0%B9_%D0%BA%D0%BE%D1%80%D0%B5%D0%BD%D1%8C_%D0%B8%D0%B7_%D0%B5%D0%B4%D0%B8%D0%BD%D0%B8%D1%86%D1%8B

[fft-online]: http://www.engineeringproductivitytools.com/stuff/T0001/index.html "Let's shed some light on the FFT"
[PFA]: https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm

[fast polynomial multiplication]: https://www.google.com/search?q=fast+polynomial+multiplication "fast polynomial multiplication"
[fast polynomial division]: https://www.google.com/search?q=fast+polynomial+division "fast polynomial division"
[fast polynomial evaluation]: https://www.google.com/search?q=fast+polynomial+evaluation "fast polynomial evaluation"
[fast polynomial interpolation]: https://www.google.com/search?q=fast+polynomial+interpolation "fast polynomial interpolation"

[Galua field]: https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D0%BD%D0%B5%D1%87%D0%BD%D0%BE%D0%B5_%D0%BF%D0%BE%D0%BB%D0%B5