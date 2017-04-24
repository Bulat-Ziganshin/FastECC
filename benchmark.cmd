@set N=%1
@if %N%.==. set N=16
@set SIZE=%2
@if %SIZE%.==. set SIZE=8192
@for %%p in ("" "=" "-" "+") do @(
    for %%e in (64g 64m 32g 32m) do @(
        echo ntt%%e %%po %N% %SIZE%
        ntt%%e %%po %N% %SIZE%
        echo.
        echo ntt%%e %%pn %N% %SIZE%
        ntt%%e %%pn %N% %SIZE%
        echo.
        echo ntt%%e %%pb
        ntt%%e %%pb
        echo.
    )
)    

@for %e% in (64g 64m 32g 32m) do rs%e 19 4000
