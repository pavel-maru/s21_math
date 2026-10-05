#include <mpfr.h>
#include <stdio.h>

int main(void) {
  mpfr_set_default_prec(256);
  mpfr_t pi, two_pi, ln2;
  mpfr_inits(pi, two_pi, ln2, (mpfr_ptr)0);
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
  mpfr_const_log2(ln2, MPFR_RNDN);

  mpfr_t hi_113, lo;
  mpfr_init2(hi_113, 113);
  mpfr_init2(lo, 256);

  mpfr_set(hi_113, two_pi, MPFR_RNDN);
  mpfr_sub(lo, two_pi, hi_113, MPFR_RNDN);
  mpfr_printf("2PI_HI = %.40Rg\n2PI_LO = %.40Rg\n", hi_113, lo);

  mpfr_set(hi_113, ln2, MPFR_RNDN);
  mpfr_sub(lo, ln2, hi_113, MPFR_RNDN);
  mpfr_printf("LN2_HI = %.40Rg\nLN2_LO = %.40Rg\n", hi_113, lo);

  mpfr_clears(pi, two_pi, ln2, hi_113, lo, (mpfr_ptr)0);
  return 0;
}
