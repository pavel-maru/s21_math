#include "s21_math.h"

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

#define S21_EPS 1e-25L
#define S21_MAX_ITER 300

#define S21_INF __builtin_infl()
#define S21_NAN __builtin_nanl("")

#define S21_LL_MAX 9223372036854775807.0L
#define S21_LL_MIN -9223372036854775808.0L

/* exp(y) с |y| > ~11356 переполняет 80-битный long double */
#define S21_EXP_LIMIT 11356.0L

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

/* x * 2^n для неотрицательного n (n >= 0 гарантируется вызывающим) */
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

/* Приведение x по модулю 2*pi с использованием разбитой константы */
static long double s21_reduce_2pi(long double x) {
  long double n = s21_trunc_l(x / S21_2PI);
  long double r = ((x - n * S21_2PI_HI) - n * S21_2PI_LO);
  return r;
}

/* Kahan summation: sum += term с компенсацией c */
#define S21_KAHAN_ADD(sum, c, term) \
  do {                              \
    long double _y = (term) - (c);  \
    long double _t = (sum) + _y;    \
    (c) = (_t - (sum)) - _y;        \
    (sum) = _t;                     \
  } while (0)

/* Бинарное возведение в степень: base^n, n - целое */
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

  /* Знак результата = знак x; работаем с модулями.
     Внутренний цикл подбирает максимальное t = y * 2^k <= xl,
     внешний вычитает такие t, пока xl >= yl.
     Это устойчиво к большим x, где xl / yl не влезает в long long. */
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

  /* Newton–Raphson: квадратичная сходимость */
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

  if (y > S21_EXP_LIMIT) {
    return sign ? 0.0L : S21_INF;
  }

  /* exp(y) = 2^k * exp(r), y = k*ln2 + r, |r| <= ln2/2 */
  long double k = s21_trunc_l(y / S21_LN2 + 0.5L);
  long double r = y - k * S21_LN2_HI - k * S21_LN2_LO;

  long double term = 1.0L;
  long double sum = 1.0L, c = 0.0L;
  for (int i = 1; i < S21_MAX_ITER; i++) {
    term *= r / (long double)i;
    S21_KAHAN_ADD(sum, c, term);
    if (s21_fabsl(term) < S21_EPS) break;
  }

  long double result = s21_ldexp_int(sum, (long long)k);
  if (sign) result = 1.0L / result;
  return result;
}

long double s21_log(double x) {
  if (s21_isnan_l(x)) return x;
  if (x < 0) return S21_NAN;
  if (x == 0) return -S21_INF;
  if (s21_isinf_l(x)) return x;
  if (x == 1.0) return 0.0L;

  long double y = (long double)x;
  long long count = 0;

  /* Приведение y к [1, 2) */
  while (y >= 2.0L) {
    y /= 2.0L;
    count++;
  }
  while (y < 1.0L) {
    y *= 2.0L;
    count--;
  }

  /* z = (y-1)/(y+1) ∈ [0, 1/3) -> быстрая сходимость */
  long double z = (y - 1.0L) / (y + 1.0L);
  long double z2 = z * z;
  long double term = z;
  long double sum = 0.0L, c = 0.0L;

  for (int i = 1; i < S21_MAX_ITER; i += 2) {
    S21_KAHAN_ADD(sum, c, term / (long double)i);
    term *= z2;
    if (s21_fabsl(term) < S21_EPS) break;
  }

  /* 2*sum + count*ln2, ln2 разбит на HI + LO */
  long double result = 2.0L * sum + (long double)count * S21_LN2_HI +
                       (long double)count * S21_LN2_LO;
  return result;
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

  /* Целый показатель — точное бинарное возведение, без exp/log.
     Порог 1e18 — защита от переполнения long long при приведении. */
  long double e_int = s21_trunc_l(e);
  if (e_int == e && s21_fabsl(e) < 1e18L) {
    return s21_powi(b, (long long)e_int);
  }

  /* Дробный показатель при отрицательном основании — NaN */
  if (b < 0.0L) return S21_NAN;

  return s21_exp(e * s21_log(b));
}

long double s21_sin(double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return S21_NAN;

  long double xl = s21_reduce_2pi((long double)x);

  long double term = xl;
  long double sum = xl, c = 0.0L;
  long double x2 = xl * xl;

  for (int i = 3; i < S21_MAX_ITER; i += 2) {
    term *= -x2 / ((long double)(i - 1) * (long double)i);
    S21_KAHAN_ADD(sum, c, term);
    if (s21_fabsl(term) < S21_EPS) break;
  }
  return sum;
}

long double s21_cos(double x) {
  if (s21_isnan_l(x) || s21_isinf_l(x)) return S21_NAN;

  long double xl = s21_reduce_2pi((long double)x);

  long double term = 1.0L;
  long double sum = 1.0L, c = 0.0L;
  long double x2 = xl * xl;

  for (int i = 2; i < S21_MAX_ITER; i += 2) {
    term *= -x2 / ((long double)(i - 1) * (long double)i);
    S21_KAHAN_ADD(sum, c, term);
    if (s21_fabsl(term) < S21_EPS) break;
  }
  return sum;
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

  /* Argument reduction: atan(x) = 2*atan(x / (1 + sqrt(1+x^2))) */
  int reductions = 0;
  while (xl > 0.1L && reductions < 8) {
    xl = xl / (1.0L + s21_sqrt(1.0L + xl * xl));
    reductions++;
  }

  long double term = xl;
  long double sum = xl, c = 0.0L;
  long double x2 = xl * xl;

  for (int i = 3; i < S21_MAX_ITER; i += 2) {
    term *= -x2;
    S21_KAHAN_ADD(sum, c, term / (long double)i);
    if (s21_fabsl(term / (long double)i) < S21_EPS) break;
  }

  for (int i = 0; i < reductions; i++) sum *= 2.0L;

  if (invert) sum = S21_PI_2 - sum;
  if (sign) sum = -sum;
  return sum;
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
