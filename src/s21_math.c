#include "s21_math.h"

/* Определяем разрядность long double.
   __LDBL_MANT_DIG__ значения:
     53  — 64-бит (Apple Silicon, ARM64 macOS, -mlong-double-64)
     64  — 80-бит (x86-64 extended precision)
     113 — 128-бит (IEEE 754 quad, ARM64 Linux) */
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
#define S21_LDBL_IS_DOUBLE 1
#else
#define S21_LDBL_IS_DOUBLE 0
#endif

#define S21_PI 3.1415926535897932384626433832795028841971693993751L
#define S21_PI_2 1.5707963267948966192313216916397514420985846996876L
#define S21_2PI 6.2831853071795864769252867665590057683943387987502L

/* 2*pi = HI + LO, чтобы точнее приводить аргумент по модулю 2*pi */
#define S21_2PI_HI 6.283185307179586L
#define S21_2PI_LO 2.4492935982947064e-16L

/* ln2 = HI + LO для exp/log */
#define S21_LN2 0.6931471805599453094172321214581765680755001343603L
#define S21_LN2_HI 0.693147180559945309417232121458176568L
#define S21_LN2_LO 7.5500134360255254121e-33L

/* Параметры зависят от точности long double.
   На Apple Silicon long double == double, эпсилон ~2.2e-16,
   поэтому порог сходимости поднят до 1e-15, а лимит exp
   снижен до ln(DBL_MAX) ≈ 709.78 (переполнение наступает раньше). */
#if S21_LDBL_IS_DOUBLE
#define S21_EPS 1e-15L
#define S21_MAX_ITER 100
#define S21_EXP_LIMIT 709.0L
#else
#define S21_EPS 1e-25L
#define S21_MAX_ITER 300
#define S21_EXP_LIMIT 11356.0L
#endif

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

static long double s21_trunc_l(long double x) {
  if (x >= 0) {
    if (x >= S21_LL_MAX) return x;
    return (long double)(long long)x;
  }
  if (x <= S21_LL_MIN) return x;
  return (long double)(long long)x;
}

/* x * 2^n для неотрицательного n */
static long double s21_ldexp_int(long double x, long long n) {
  if (n == 0) return x;
  long double result = x;
  long double factor = 2.0L;
  while (n) {
    if (n & 1) result *= factor;
    factor *= factor;
    n >>= 1;
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

/* Точное сложение a + b: s + err == a + b (по математике). */
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

/* Точное произведение a·b: p + err == a·b (по математике).
   Использует FMA / fma. В 64-битном режиме (Apple Silicon)
   использует fma для double, так как fmal несовместима с ABI
   glibc при -mlong-double-64. */
static void s21_two_product(long double a, long double b, long double *p,
                            long double *err) {
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  /* long double == double: fma для double */
  double a_d = (double)a;
  double b_d = (double)b;
  double p_d = a_d * b_d;
  double err_d = __builtin_fma(a_d, b_d, -p_d);
  *p = (long double)p_d;
  *err = (long double)err_d;
#else
  /* 80-бит или 128-бит: fmal работает корректно */
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

/* a + b с сохранением потерянных битов и ре-нормализацией. */
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

/* Приведение x по модулю 2π с использованием double-double.

   Ключевая идея: вычитание x − n·2π выполняется в dd, чтобы
   сохранить младшие биты. Это критично для sin/cos на больших
   аргументах, где x и n·2π близки, а разность мала.

   Раньше вычитание шло в long double, что давало ошибку
   ~|x|·ε ≈ 1e10·1e-16 = 1e-6 в остатке r и приводило к ULP ~1e9
   на sin(1e10) в 64-битном режиме. Теперь ошибка ~1e10·ε² ≈ 1e-22. */
static long double s21_reduce_2pi(long double x) {
  /* n = trunc(x / 2π) */
  long double n = s21_trunc_l(x / S21_2PI);

  /* Вычислим n·2π в dd: n·2π = n·2π_hi + n·2π_lo
     n·2π_hi — точное произведение через two_product */
  long double p_hi, p_err;
  s21_two_product(n, S21_2PI_HI, &p_hi, &p_err);
  long double p_lo = n * S21_2PI_LO;

  s21_dd n_2pi;
  s21_two_sum(p_hi, p_err + p_lo, &n_2pi.hi, &n_2pi.lo);

  /* r = x − n·2π в dd */
  s21_dd x_dd = s21_dd_make(x);
  s21_dd neg;
  neg.hi = -n_2pi.hi;
  neg.lo = -n_2pi.lo;
  s21_dd r_dd = s21_dd_add(x_dd, neg);

  return s21_dd_value(r_dd);
}

/* ============================================================
   Бинарное возведение в степень
   ============================================================ */

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
   Публичные функции
   ============================================================ */

int s21_abs(int x) { return x < 0 ? -x : x; }

long double s21_fabs(double x) {
  if (x < 0) return -(long double)x;
  if (x == 0) return 0.0L;
  return (long double)x;
}

long double s21_ceil(double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return x;
  long double y = (long double)x;
  if (y >= S21_LL_MAX || y <= S21_LL_MIN) return y;
  long long i = (long long)y;
  if (y > 0 && y > (long double)i) return (long double)(i + 1);
  if (y < 0 && y < (long double)i) return (long double)i;
  return (long double)i;
}

long double s21_floor(double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return x;
  long double y = (long double)x;
  if (y >= S21_LL_MAX || y <= S21_LL_MIN) return y;
  long long i = (long long)y;
  if (y > 0 && y > (long double)i) return (long double)i;
  if (y < 0 && y < (long double)i) return (long double)(i - 1);
  return (long double)i;
}

long double s21_fmod(double x, double y) {
  if (s21_isnan_l(x) || s21_isnan_l(y) || s21_isinf_l(x) || y == 0.0)
    return S21_NAN;
  if (s21_isinf_l(y)) return x;

  long double xl = s21_fabsl((long double)x);
  long double yl = s21_fabsl((long double)y);

  while (xl >= yl) {
    long double t = yl;
    while (t * 2.0L <= xl) t *= 2.0L;
    xl -= t;
  }

  return ((long double)x < 0) ? -xl : xl;
}

long double s21_sqrt(double x) {
  if (s21_isnan_l(x)) return x;
  if (x < 0) return S21_NAN;
  if (x == 0 || s21_isinf_l(x)) return x;

  long double xl = (long double)x;
  long double res = (xl < 1.0L) ? 1.0L : xl;

  for (int i = 0; i < 100; i++) {
    long double next = 0.5L * (res + xl / res);
    if (next == res) break;
    res = next;
  }
  return res;
}

long double s21_exp(double x) {
  if (s21_isnan_l(x)) return x;
  if (s21_isinf_l(x)) return x > 0 ? S21_INF : 0.0L;
  if (x == 0.0) return 1.0L;

  long double y = (long double)x;
  int sign = 0;
  if (y < 0) {
    sign = 1;
    y = -y;
  }
  if (y > S21_EXP_LIMIT) return sign ? 0.0L : S21_INF;

  long double k = s21_trunc_l(y / S21_LN2 + 0.5L);

  long double r;

#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 64
  /* 80-бит: FMA-путь через __builtin_fmal */
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
  /* 128-бит или 64-бит: обычное разбитое вычитание */
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

long double s21_log(double x) {
  if (s21_isnan_l(x)) return x;
  if (x < 0) return S21_NAN;
  if (x == 0) return -S21_INF;
  if (s21_isinf_l(x)) return x;
  if (x == 1.0) return 0.0L;

  long double y = (long double)x;
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

  s21_dd t1 = s21_dd_make((long double)count * S21_LN2_HI);
  s21_dd t2 = s21_dd_make((long double)count * S21_LN2_LO);

  s21_dd result = s21_dd_add(twice, t1);
  result = s21_dd_add(result, t2);
  return s21_dd_value(result);
}

long double s21_pow(double base, double exp_val) {
  if (s21_isnan_l(base) || s21_isnan_l(exp_val)) return S21_NAN;
  if (exp_val == 0.0) return 1.0L;

  long double b = (long double)base;
  long double e = (long double)exp_val;

  if (b == 0.0L) {
    if (e < 0.0L) return S21_INF;
    return 0.0L;
  }

  long double e_int = s21_trunc_l(e);
  if (e_int == e && s21_fabsl(e) < 1e18L) {
    return s21_powi(b, (long long)e_int);
  }

  if (b < 0.0L) return S21_NAN;

  return s21_exp(e * s21_log(b));
}

long double s21_sin(double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return S21_NAN;

  long double xl = s21_reduce_2pi((long double)x);

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

long double s21_cos(double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return S21_NAN;

  long double xl = s21_reduce_2pi((long double)x);

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

long double s21_tan(double x) {
  long double s = s21_sin(x);
  long double c = s21_cos(x);
  if (c == 0.0L) return S21_INF;
  return s / c;
}

long double s21_atan(double x) {
  if (s21_isnan_l(x)) return x;
  if (s21_isinf_l(x)) return x > 0 ? S21_PI_2 : -S21_PI_2;

  long double xl = (long double)x;
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
    xl = xl / (1.0L + s21_sqrt(1.0L + xl * xl));
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

long double s21_asin(double x) {
  if (s21_isnan_l(x)) return x;
  if (x < -1.0 || x > 1.0) return S21_NAN;
  if (x == 1.0) return S21_PI_2;
  if (x == -1.0) return -S21_PI_2;
  return s21_atan((long double)x / s21_sqrt(1.0 - (long double)x * x));
}

long double s21_acos(double x) {
  if (s21_isnan_l(x)) return x;
  if (x < -1.0 || x > 1.0) return S21_NAN;
  return S21_PI_2 - s21_asin(x);
}
