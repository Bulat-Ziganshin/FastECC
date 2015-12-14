#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>

typedef unsigned long long NUM;

bool is_prime(NUM N)
{
  NUM m = std::min(NUM(sqrt(N))+10, N);
  for (NUM i=3; i < m; i+=2)
    if (N%2==0 || N%i==0)
      return false;
  return true;
}

NUM max_divider (NUM N)
{
  NUM m = std::min(NUM(sqrt(N))+10, N), N0 = N, max_divider = 0;
  std::cout << N << " = ";
  for (NUM i=2; N>1 && i<m; i++)
    for ( ; N%i==0; N/=i, max_divider = i)
      std::cout << (max_divider? "*":"") << i;
  if (N==1) {
    return max_divider;
  } else if (N<N0) {
    std::cout << (max_divider? "*":"") << N;
    return N;
  } else {
    return N;
  }
}


int main (int argc, char **argv)
{
  if(argc!=2)  {printf("Usage: prime N\n  Prints prime number higher or equal to N\n"); return 0;}
  NUM N, N0;
  sscanf (argv[1], "%llu", &N);  N0=N;
  if (max_divider(N) < N)
    for (printf("\n"); N++; )
    {
      NUM m = std::min(NUM(sqrt(N))+10, N);
      for (NUM i=3; i < m; i+=2)
        if (N%2==0 || N%i==0)
          {if (N%2==0)  i=2;  printf("%llu / %u%s", N, i, (i>1000 || N==N0)? "\n" : "\r"); goto next;}
      break;
      next: ;
    }
  std::cout << "\r" << N << " is prime!!!\n";

  NUM max_d = N;
  for (N=N0; max_d>16; N++)
  {
    NUM d = max_divider(N);
    if (d!=N  &&  d<max_d  &&  is_prime(N+1)) {
      max_d = d;
      std::cout << "\n";
    } else {
      std::cout << "\r                                                                  \r";
    }
  }

  return 0;
}
