/* probe_consts2.c — исправленный probe констант.
   Печатает HI и LO для 2π и ln2 в текущем режиме long double.
   Использует __float128 для точного вычисления остатка. */

#include <stdio.h>
#include <math.h>

#if defined(__SIZEOF_FLOAT128__)
typedef __float128 quad_t;
#define Q(x) x##Q
#else
typedef long double quad_t;
#define Q(x) x##L
#endif

static void print_ld(long double x, const char *tag) {
  printf("%s = ", tag);
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  /* 64-бит: %Le сломан, печатаем как double */
  printf("%.20e (double)", (double)x);
#elif defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 113
  /* 128-бит: %Le тоже сломан, печатаем hex-дамп 16 байт MSB-first */
  unsigned char *p = (unsigned char *)&x;
  for (int i = (int)sizeof(x) - 1; i >= 0; i--) printf("%02x", p[i]);
  printf(" (hex)");
#else
  /* 80-бит: %Le работает */
  printf("%.25Le", x);
#endif
  printf("\n");
}

int main(void) {
  printf("MANT_DIG=%d  sizeof(long double)=%zu  sizeof(quad_t)=%zu\n",
         __LDBL_MANT_DIG__, sizeof(long double), sizeof(quad_t));

  /* Полные значения с большим числом знаков — компилятор сам
     округлит их до текущей точности long double. */
  quad_t pi_true =
      Q(6.283185307179586476925286766559005768394338798750211641949);
  quad_t ln2_true =
      Q(0.69314718055994530941723212145817656807550013436025525412);

  long double pi_hi  = (long double)pi_true;
  long double ln2_hi = (long double)ln2_true;

  /* Точный остаток true − hi, вычисленный в __float128. */
  quad_t pi_lo_q  = pi_true  - (quad_t)pi_hi;
  quad_t ln2_lo_q = ln2_true - (quad_t)ln2_hi;

  long double pi_lo  = (long double)pi_lo_q;
  long double ln2_lo = (long double)ln2_lo_q;

  print_ld(pi_hi,  "2PI_HI");
  print_ld(pi_lo,  "2PI_LO");
  print_ld(ln2_hi, "LN2_HI");
  print_ld(ln2_lo, "LN2_LO");
  return 0;
}
