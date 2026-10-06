#include "s21_math.h"

/* Разрядность long double определяется через __LDBL_MANT_DIG__:
     53  — double (Apple Silicon, ARM64 macOS, -mlong-double-64)
     64  — 80-бит x86 extended precision
     113 — 128-бит IEEE 754 quad (ARM64 Linux, -mlong-double-128)

   Константы 2π и ln2 заданы парами HI + LO, где HI — ближайшее
   представимое long double, а LO — точный остаток (константа − HI).
   Значения получены эмпирически:
     80-бит  — probe_consts2.c
     64-бит  — probe_consts2.c
     128-бит — probe_128_mpfr.c (MPFR с 256-битной точностью) */
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  #define S21_2PI_HI  6.28318530717958623200L
  #define S21_2PI_LO  2.44929359829470641435e-16L
  #define S21_LN2_HI  0.693147180559945286227L
  #define S21_LN2_LO  2.31904681384629955842e-17L
  #define S21_EPS         1e-15L
  #define S21_MAX_ITER    100
  #define S21_EXP_LIMIT   709.0L
#elif defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 64
  #define S21_2PI_HI  6.2831853071795864770256179L
  #define S21_2PI_LO -1.0033115225336681390734914e-19L
  #define S21_LN2_HI  0.69314718055994530942869047L
  #define S21_LN2_LO -1.1458352726798725802795825e-20L
  #define S21_EPS         1e-25L
  #define S21_MAX_ITER    300
  #define S21_EXP_LIMIT   11356.0L
#elif defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 113
  #define S21_2PI_HI  6.283185307179586476925286766559005594958L
  #define S21_2PI_LO  1.734362026024756204959408805208670393752e-34L
  #define S21_LN2_HI  0.6931471805599453094172321214581765750836L
  #define S21_LN2_LO -7.008139474549585163412662008771625673778e-36L
  #define S21_EPS         1e-35L
  #define S21_MAX_ITER    400
  #define S21_EXP_LIMIT   11356.0L
#else
  #error "Unsupported __LDBL_MANT_DIG__"
#endif

/* S21_PI и S21_PI_2 — для asin/acos/atan, не участвуют в range
   reduction. S21_2PI используется только в fallback-ветке
   s21_reduce_2pi для |x| > 1.16e20. */
#define S21_PI   3.1415926535897932384626433832795028841971693993751L
#define S21_PI_2 1.5707963267948966192313216916397514420985846996876L
#define S21_2PI  6.2831853071795864769252867665590057683943387987502L

#define S21_INF __builtin_infl()
#define S21_NAN __builtin_nanl("")

#define S21_LL_MAX 9223372036854775807.0L
#define S21_LL_MIN -9223372036854775808.0L

/* ============================================================
   Базовые хелперы
   ============================================================ */

static int s21_isnan_l(long double x) { return x != x; }
static int s21_isinf_l(long double x) { return x == S21_INF || x == -S21_INF; }
static long double s21_fabsl(long double x) { return x < 0 ? -x : x; }
static int s21_signbit_l(long double x) { return __builtin_signbitl(x); }

/* Отбрасывание дробной части. Работает до 2^64 через unsigned long
   long, что расширяет корректный диапазон s21_reduce_2pi до
   |x| < 2π·2^64 ≈ 1.16e20. */
static long double s21_trunc_l(long double x) {
  if (x >= 0) {
    if (x < 9223372036854775808.0L)     /* < 2^63 */
      return (long double)(long long)x;
    if (x < 18446744073709551616.0L)    /* < 2^64 */
      return (long double)(unsigned long long)x;
    return x;
  }
  if (x > -9223372036854775808.0L)
    return (long double)(long long)x;
  if (x > -18446744073709551616.0L)
    return -(long double)(unsigned long long)(-x);
  return x;
}

/* Нечётное целое? Используется в s21_pow_l для случая base = -0.0.
   Для |e| >= 2^63 знак показателя определить нельзя без потери
   точности (long double там уже не различает чётные/нечётные),
   поэтому возвращаем 0. */
static int s21_is_odd_int(long double e) {
  long double t = s21_trunc_l(e);
  if (t != e) return 0;
  if (s21_fabsl(t) >= 9223372036854775808.0L) return 0;  /* >= 2^63 */
  return (int)(((long long)t) & 1LL);
}

/* x · 2^n для любого целого n (включая отрицательные). */
static long double s21_ldexp_int(long double x, long long n) {
  int neg = (n < 0);
  unsigned long long un =
      neg ? (unsigned long long)(-(n + 1)) + 1ULL : (unsigned long long)n;
  long double result = x;
  long double factor = neg ? 0.5L : 2.0L;
  while (un) {
    if (un & 1ULL) result *= factor;
    un >>= 1;
    if (un) factor *= factor;
  }
  return result;
}

/* ============================================================
   Расширение из 2 компонент (Shewchuk-style double-double).
   ============================================================ */

typedef struct {
  long double hi;
  long double lo;
} s21_dd;

static void s21_two_sum(long double a, long double b, long double *s,
                        long double *err) {
  long double sum = a + b;
  long double bv = sum - a;
  long double av = sum - bv;
  long double br = b - bv;
  long double ar = a - av;
  *s = sum;
  *err = ar + br;
}

static void s21_two_product(long double a, long double b, long double *p,
                            long double *err) {
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  double a_d = (double)a;
  double b_d = (double)b;
  double p_d = a_d * b_d;
  double err_d = __builtin_fma(a_d, b_d, -p_d);
  *p = (long double)p_d;
  *err = (long double)err_d;
#else
  *p = a * b;
  *err = __builtin_fmal(a, b, -*p);
#endif
}

static s21_dd s21_dd_make(long double x) {
  s21_dd r;
  r.hi = x;
  r.lo = 0.0L;
  return r;
}

static s21_dd s21_dd_add(s21_dd a, s21_dd b) {
  long double s, err;
  s21_two_sum(a.hi, b.hi, &s, &err);
  err += a.lo + b.lo;

  long double s2, err2;
  s21_two_sum(s, err, &s2, &err2);

  s21_dd r;
  r.hi = s2;
  r.lo = err2;
  return r;
}

static s21_dd s21_dd_ldexp(s21_dd a, long long n) {
  a.hi = s21_ldexp_int(a.hi, n);
  a.lo = s21_ldexp_int(a.lo, n);
  return a;
}

static long double s21_dd_value(s21_dd a) { return a.hi + a.lo; }

static long double s21_reduce_2pi(long double x) {
  if (s21_fabsl(x) > 1.16e20L) {
    long double r = x;
    while (s21_fabsl(r) >= S21_2PI) {
      long double t = S21_2PI;
      while (s21_fabsl(r) >= t * 2.0L) t *= 2.0L;
      long double prev = r;
      r = (r >= 0) ? r - t : r + t;
      if (r == prev) break;
    }
    return r;
  }

  long double n = s21_trunc_l(x / S21_2PI);

  long double p_hi, p_err;
  s21_two_product(n, S21_2PI_HI, &p_hi, &p_err);
  long double p_lo = n * S21_2PI_LO;

  s21_dd n_2pi;
  s21_two_sum(p_hi, p_err + p_lo, &n_2pi.hi, &n_2pi.lo);

  s21_dd x_dd = s21_dd_make(x);
  s21_dd neg;
  neg.hi = -n_2pi.hi;
  neg.lo = -n_2pi.lo;
  s21_dd r_dd = s21_dd_add(x_dd, neg);
  return s21_dd_value(r_dd);
}

static long double s21_powi(long double base, long long n) {
  int neg = (n < 0);
  unsigned long long m =
      neg ? (unsigned long long)(-(n + 1)) + 1ULL : (unsigned long long)n;

  long double result = 1.0L;
  long double factor = base;
  while (m) {
    if (m & 1ULL) result *= factor;
    m >>= 1;
    if (m) factor *= factor;
  }
  return neg ? 1.0L / result : result;
}

/* ============================================================
   Внутренние _l-реализации.

   Публичный API принимает double (по требованию задания), но
   внутри всё считается в long double. Иначе при вызове
   s21_sqrt_l(1.0L + xl*xl) из atan аргумент усекается до double,
   и ошибка растёт на много порядков.

   Публичные функции ниже — тонкие обёртки над этими _l-версиями.
   ============================================================ */

static long double s21_fabs_l(long double x) {
  if (x < 0) return -x;
  if (x == 0) return 0.0L;
  return x;
}

static long double s21_ceil_l(long double y) {
  if (s21_isnan_l(y) || s21_isinf_l(y)) return y;
  if (y == 0.0L) return y;
  if (y >= S21_LL_MAX || y <= S21_LL_MIN) return y;
  long long i = (long long)y;
  if (y > 0 && y > (long double)i) return (long double)(i + 1);
  if (y < 0 && y < (long double)i) return (long double)i;
  return (long double)i;
}

static long double s21_floor_l(long double y) {
  if (s21_isnan_l(y) || s21_isinf_l(y)) return y;
  if (y == 0.0L) return y;
  if (y >= S21_LL_MAX || y <= S21_LL_MIN) return y;
  long long i = (long long)y;
  if (y > 0 && y > (long double)i) return (long double)i;
  if (y < 0 && y < (long double)i) return (long double)(i - 1);
  return (long double)i;
}

static long double s21_fmod_l(long double x, long double y) {
  if (s21_isnan_l(x) || s21_isnan_l(y) || s21_isinf_l(x) || y == 0.0L)
    return S21_NAN;
  if (s21_isinf_l(y)) return x;
  if (x == 0.0L) return x;

  long double xl = s21_fabsl(x);
  long double yl = s21_fabsl(y);

  while (xl >= yl) {
    long double t = yl;
    while (t * 2.0L <= xl) t *= 2.0L;
    xl -= t;
  }

  return (x < 0) ? -xl : xl;
}

static long double s21_sqrt_l(long double x) {
  if (s21_isnan_l(x)) return x;
  if (x < 0) return S21_NAN;
  if (x == 0 || s21_isinf_l(x)) return x;

  int scale = 0;
  while (x < 0.25L) { x *= 4.0L; scale++; }
  while (x > 1.0L)  { x *= 0.25L; scale--; }

  long double res = 0.75L;
  for (int i = 0; i < 60; i++) {
    long double next = 0.5L * (res + x / res);
    if (next == res) break;
    res = next;
  }

  return s21_ldexp_int(res, -(long long)scale);
}

static long double s21_exp_l(long double x) {
  if (s21_isnan_l(x)) return x;
  if (s21_isinf_l(x)) return x > 0 ? S21_INF : 0.0L;
  if (x == 0.0L) return 1.0L;

  long double y = x;
  int sign = 0;
  if (y < 0) {
    sign = 1;
    y = -y;
  }
  if (y > S21_EXP_LIMIT) return sign ? 0.0L : S21_INF;

  long double k = s21_trunc_l(y / S21_LN2_HI + 0.5L);

  long double r;

#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 64
  long double khi = (long double)k * S21_LN2_HI;
  long double khi_err = __builtin_fmal((long double)k, S21_LN2_HI, -khi);

  long double klo = (long double)k * S21_LN2_LO;
  long double klo_err = __builtin_fmal((long double)k, S21_LN2_LO, -klo);

  s21_dd r_dd = s21_dd_make(y);
  s21_dd neg1;
  neg1.hi = -khi;
  neg1.lo = -khi_err;
  s21_dd neg2;
  neg2.hi = -klo;
  neg2.lo = -klo_err;
  r_dd = s21_dd_add(r_dd, neg1);
  r_dd = s21_dd_add(r_dd, neg2);
  r = s21_dd_value(r_dd);
#else
  r = y - k * S21_LN2_HI - k * S21_LN2_LO;
#endif

  long double term = 1.0L;
  s21_dd sum = s21_dd_make(1.0L);

  for (int i = 1; i < S21_MAX_ITER; i++) {
    term *= r / (long double)i;
    sum = s21_dd_add(sum, s21_dd_make(term));
    if (s21_fabsl(term) < S21_EPS) break;
  }

  s21_dd res = s21_dd_ldexp(sum, (long long)k);
  long double out = s21_dd_value(res);
  if (sign) out = 1.0L / out;
  return out;
}

static long double s21_log_l(long double x) {
  if (s21_isnan_l(x)) return x;
  if (x < 0) return S21_NAN;
  if (x == 0) return -S21_INF;
  if (s21_isinf_l(x)) return x;
  if (x == 1.0L) return 0.0L;

  long double y = x;
  long long count = 0;

  while (y >= 2.0L) {
    y /= 2.0L;
    count++;
  }
  while (y < 1.0L) {
    y *= 2.0L;
    count--;
  }

  long double z = (y - 1.0L) / (y + 1.0L);
  long double z2 = z * z;
  long double term = z;
  s21_dd sum = s21_dd_make(0.0L);

  for (int i = 1; i < S21_MAX_ITER; i += 2) {
    sum = s21_dd_add(sum, s21_dd_make(term / (long double)i));
    term *= z2;
    if (s21_fabsl(term) < S21_EPS) break;
  }

  s21_dd twice = s21_dd_make(2.0L * sum.hi);
  twice.lo = 2.0L * sum.lo;

  long double khi, khi_err;
  s21_two_product((long double)count, S21_LN2_HI, &khi, &khi_err);
  long double klo = (long double)count * S21_LN2_LO;

  s21_dd result = twice;
  result = s21_dd_add(result, s21_dd_make(khi));
  result = s21_dd_add(result, s21_dd_make(khi_err));
  result = s21_dd_add(result, s21_dd_make(klo));
  return s21_dd_value(result);
}

static long double s21_pow_l(long double b, long double e) {
  if (s21_isnan_l(b) || s21_isnan_l(e)) return S21_NAN;
  if (e == 0.0L) return 1.0L;

  if (b == 0.0L) {
    int b_neg = s21_signbit_l(b);
    int e_odd = s21_is_odd_int(e);

    if (e < 0.0L) return (b_neg && e_odd) ? -S21_INF : S21_INF;
    return (b_neg && e_odd) ? -0.0L : 0.0L;
  }

  long double e_int = s21_trunc_l(e);
  if (e_int == e && s21_fabsl(e) < 1e18L) {
    return s21_powi(b, (long long)e_int);
  }

  if (b < 0.0L) return S21_NAN;

  return s21_exp_l(e * s21_log_l(b));
}

static long double s21_sin_l(long double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return S21_NAN;

  long double xl = s21_reduce_2pi(x);

  long double term = xl;
  s21_dd sum = s21_dd_make(xl);
  long double x2 = xl * xl;

  for (int i = 3; i < S21_MAX_ITER; i += 2) {
    term *= -x2 / ((long double)(i - 1) * (long double)i);
    sum = s21_dd_add(sum, s21_dd_make(term));
    if (s21_fabsl(term) < S21_EPS) break;
  }
  return s21_dd_value(sum);
}

static long double s21_cos_l(long double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return S21_NAN;

  long double xl = s21_reduce_2pi(x);

  long double term = 1.0L;
  s21_dd sum = s21_dd_make(1.0L);
  long double x2 = xl * xl;

  for (int i = 2; i < S21_MAX_ITER; i += 2) {
    term *= -x2 / ((long double)(i - 1) * (long double)i);
    sum = s21_dd_add(sum, s21_dd_make(term));
    if (s21_fabsl(term) < S21_EPS) break;
  }
  return s21_dd_value(sum);
}

static long double s21_tan_l(long double x) {
  long double s = s21_sin_l(x);
  long double c = s21_cos_l(x);
  if (c == 0.0L) return S21_INF;
  return s / c;
}

static long double s21_atan_l(long double xl) {
  if (s21_isnan_l(xl)) return xl;
  if (s21_isinf_l(xl)) return xl > 0 ? S21_PI_2 : -S21_PI_2;

  int sign = 0;
  if (xl < 0) {
    sign = 1;
    xl = -xl;
  }

  int invert = 0;
  if (xl > 1.0L) {
    invert = 1;
    xl = 1.0L / xl;
  }

  int reductions = 0;
  while (xl > 0.1L && reductions < 8) {
    xl = xl / (1.0L + s21_sqrt_l(1.0L + xl * xl));
    reductions++;
  }

  long double term = xl;
  s21_dd sum = s21_dd_make(xl);
  long double x2 = xl * xl;

  for (int i = 3; i < S21_MAX_ITER; i += 2) {
    term *= -x2;
    sum = s21_dd_add(sum, s21_dd_make(term / (long double)i));
    if (s21_fabsl(term / (long double)i) < S21_EPS) break;
  }

  for (int i = 0; i < reductions; i++) {
    sum.hi *= 2.0L;
    sum.lo *= 2.0L;
  }

  if (invert) {
    s21_dd pi2 = s21_dd_make(S21_PI_2);
    s21_dd neg;
    neg.hi = -sum.hi;
    neg.lo = -sum.lo;
    sum = s21_dd_add(pi2, neg);
  }

  long double out = s21_dd_value(sum);
  return sign ? -out : out;
}

static long double s21_asin_l(long double x) {
  if (s21_isnan_l(x)) return x;
  if (x < -1.0L || x > 1.0L) return S21_NAN;
  if (x == 1.0L) return S21_PI_2;
  if (x == -1.0L) return -S21_PI_2;
  return s21_atan_l(x / s21_sqrt_l(1.0L - x * x));
}

static long double s21_acos_l(long double x) {
  if (s21_isnan_l(x)) return x;
  if (x < -1.0L || x > 1.0L) return S21_NAN;
  return S21_PI_2 - s21_asin_l(x);
}

/* ============================================================
   Публичные функции — тонкие обёртки над _l-версиями.
   Параметры double сохранены по требованию School 21.
   ============================================================ */

int s21_abs(int x) {
  if (x >= 0) return x;
  return (int)(-(unsigned)x);
}

long double s21_fabs(double x) { return s21_fabs_l((long double)x); }
long double s21_ceil(double x) { return s21_ceil_l((long double)x); }
long double s21_floor(double x) { return s21_floor_l((long double)x); }
long double s21_fmod(double x, double y) {
  return s21_fmod_l((long double)x, (long double)y);
}
long double s21_sqrt(double x) { return s21_sqrt_l((long double)x); }
long double s21_exp(double x) { return s21_exp_l((long double)x); }
long double s21_log(double x) { return s21_log_l((long double)x); }
long double s21_pow(double base, double exp_val) {
  return s21_pow_l((long double)base, (long double)exp_val);
}
long double s21_sin(double x) { return s21_sin_l((long double)x); }
long double s21_cos(double x) { return s21_cos_l((long double)x); }
long double s21_tan(double x) { return s21_tan_l((long double)x); }
long double s21_atan(double x) { return s21_atan_l((long double)x); }
long double s21_asin(double x) { return s21_asin_l((long double)x); }
long double s21_acos(double x) { return s21_acos_l((long double)x); }
