#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_math.h"

#define S21_PI 3.14159265358979323846
#define S21_PI_2 1.57079632679489661923

/* Допуск зависит от платформы:
     53  — 64-бит long double (Apple Silicon, ~2.2e-16)
     64  — 80-бит long double (x86-64, ~1.1e-19)
     113 — 128-бит long double (ARM64 Linux, ~1.9e-34) */
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
#define EPS 1e-10
#else
#define EPS 1e-12
#endif

#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
#define S21_NEXTAFTER(x, y) nextafter((double)(x), (double)(y))
#define S21_LIBM_FN(name, x) ((long double)name((double)(x)))
#else
#define S21_NEXTAFTER(x, y) nextafterl((x), (y))
#define S21_LIBM_FN(name, x) name##l(x)
#endif

#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
#define S21_PRINT_LD(x) printf("%-11.2e", (double)(x))
#else
#define S21_PRINT_LD(x) printf("%-11.2Le", (x))
#endif

/* Абсолютное отклонение. Если обе стороны равны (в том числе обе
   INF или обе NaN), считаем отклонение нулевым — иначе |INF−INF|
   даёт nan в таблице. */
#define S21_ABS_DIFF(a, b) (((a) == (b)) ? 0.0L : fabsl((a) - (b)))

/* ---------- abs ---------- */
START_TEST(test_abs_basic) {
  ck_assert_int_eq(s21_abs(0), 0);
  ck_assert_int_eq(s21_abs(5), 5);
  ck_assert_int_eq(s21_abs(-5), 5);
  ck_assert_int_eq(s21_abs(2147483647), 2147483647);
  /* INT_MIN: значение не определено стандартом, но не должно быть UB
     и должно совпадать с поведением glibc (INT_MIN на two's complement). */
  ck_assert_int_eq(s21_abs(-2147483647 - 1), -2147483647 - 1);
}
END_TEST

/* ---------- fabs ---------- */
START_TEST(test_fabs_loop) {
  double x = _i * 0.5 - 10.0;
  ck_assert_ldouble_eq_tol(s21_fabs(x), fabs(x), EPS);
}
END_TEST

/* ---------- ceil ---------- */
START_TEST(test_ceil_loop) {
  double x = _i * 0.37 - 15.0;
  ck_assert_ldouble_eq_tol(s21_ceil(x), ceil(x), EPS);
}
END_TEST

START_TEST(test_ceil_edge) {
  ck_assert_ldouble_eq(s21_ceil(0.0), ceil(0.0));
  ck_assert_ldouble_eq(s21_ceil(-0.0), ceil(-0.0));
  ck_assert_ldouble_eq(s21_ceil(1.0), ceil(1.0));
  ck_assert_ldouble_eq(s21_ceil(-1.0), ceil(-1.0));
  ck_assert(s21_ceil(INFINITY) == INFINITY);
  ck_assert(s21_ceil(-INFINITY) == -INFINITY);
}
END_TEST

/* ---------- floor ---------- */
START_TEST(test_floor_loop) {
  double x = _i * 0.41 - 15.0;
  ck_assert_ldouble_eq_tol(s21_floor(x), floor(x), EPS);
}
END_TEST

START_TEST(test_floor_edge) {
  ck_assert_ldouble_eq(s21_floor(0.0), floor(0.0));
  ck_assert_ldouble_eq(s21_floor(-0.0), floor(-0.0));
  ck_assert_ldouble_eq(s21_floor(1.0), floor(1.0));
  ck_assert_ldouble_eq(s21_floor(-1.0), floor(-1.0));
  ck_assert(s21_floor(INFINITY) == INFINITY);
  ck_assert(s21_floor(-INFINITY) == -INFINITY);
}
END_TEST

/* ---------- fmod ---------- */
START_TEST(test_fmod_loop) {
  double x = _i * 1.7 - 20.0;
  double y = 3.3;
  ck_assert_ldouble_eq_tol(s21_fmod(x, y), fmod(x, y), EPS);
}
END_TEST

START_TEST(test_fmod_edge) {
  ck_assert(s21_fmod(1.0, 0.0) != s21_fmod(1.0, 0.0));
  ck_assert_ldouble_eq(s21_fmod(5.0, INFINITY), 5.0);
  ck_assert(s21_fmod(INFINITY, 2.0) != s21_fmod(INFINITY, 2.0));
  /* fmod(±0.0, y) сохраняет знак нуля. */
  ck_assert_ldouble_eq(s21_fmod(0.0, 3.0), fmod(0.0, 3.0));
  ck_assert_ldouble_eq(s21_fmod(-0.0, 3.0), fmod(-0.0, 3.0));
}
END_TEST

START_TEST(test_fmod_large) {
  ck_assert_ldouble_eq_tol(s21_fmod(1e18, 7.0), fmod(1e18, 7.0), EPS);
  ck_assert_ldouble_eq_tol(s21_fmod(-1e18, 7.0), fmod(-1e18, 7.0), EPS);
  ck_assert_ldouble_eq_tol(s21_fmod(1e30, 3.0), fmod(1e30, 3.0), 1e-6);
  ck_assert_ldouble_eq_tol(s21_fmod(-1e30, 3.0), fmod(-1e30, 3.0), 1e-6);
}
END_TEST

/* ---------- sqrt ---------- */
START_TEST(test_sqrt_loop) {
  double x = _i * 0.7;
  ck_assert_ldouble_eq_tol(s21_sqrt(x), sqrt(x), EPS);
}
END_TEST

START_TEST(test_sqrt_edge) {
  ck_assert_ldouble_eq(s21_sqrt(0.0), 0.0);
  ck_assert_ldouble_eq(s21_sqrt(1.0), 1.0);
  ck_assert(s21_sqrt(-1.0) != s21_sqrt(-1.0));
  ck_assert(s21_sqrt(INFINITY) == INFINITY);
  /* Очень малые x: масштабирование в sqrt должно давать корректный
     результат за 60 итераций. */
  ck_assert_ldouble_eq_tol(s21_sqrt(1e-30), sqrt(1e-30), 1e-22);
  ck_assert_ldouble_eq_tol(s21_sqrt(1e-300), sqrt(1e-300), 1e-160);
}
END_TEST

/* ---------- exp ---------- */
START_TEST(test_exp_loop) {
  double x = _i * 0.3 - 6.0;
  ck_assert_ldouble_eq_tol(s21_exp(x), exp(x), EPS);
}
END_TEST

START_TEST(test_exp_edge) {
  ck_assert_ldouble_eq(s21_exp(0.0), 1.0);
  ck_assert(s21_exp(INFINITY) == INFINITY);
  ck_assert(s21_exp(-INFINITY) == 0.0);
  ck_assert(s21_exp(1e5) == INFINITY);
  ck_assert(s21_exp(-1e5) == 0.0);
}
END_TEST

/* ---------- log ---------- */
START_TEST(test_log_loop) {
  double x = _i * 0.5 + 0.5;
  ck_assert_ldouble_eq_tol(s21_log(x), log(x), EPS);
}
END_TEST

START_TEST(test_log_edge) {
  ck_assert_ldouble_eq(s21_log(1.0), 0.0);
  ck_assert(s21_log(0.0) == -INFINITY);
  ck_assert(s21_log(-1.0) != s21_log(-1.0));
  ck_assert(s21_log(INFINITY) == INFINITY);
}
END_TEST

/* ---------- pow ---------- */
START_TEST(test_pow_loop) {
  double base = _i * 0.3 + 0.5;
  ck_assert_ldouble_eq_tol(s21_pow(base, 2.0), pow(base, 2.0), EPS);
  ck_assert_ldouble_eq_tol(s21_pow(base, 0.5), pow(base, 0.5), EPS);
}
END_TEST

START_TEST(test_pow_edge) {
  ck_assert_ldouble_eq(s21_pow(2.0, 0.0), 1.0);
  ck_assert_ldouble_eq(s21_pow(0.0, 0.0), 1.0);
  ck_assert(s21_pow(0.0, -1.0) == INFINITY);
  ck_assert_ldouble_eq(s21_pow(0.0, 1.0), 0.0);
  ck_assert_ldouble_eq(s21_pow(0.0, 5.0), 0.0);
  ck_assert_ldouble_eq_tol(s21_pow(-2.0, 3.0), -8.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(-2.0, 2.0), 4.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(2.0, 10.0), 1024.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(2.0, -2.0), 0.25, EPS);
  ck_assert(s21_pow(-2.0, 0.5) != s21_pow(-2.0, 0.5));
  /* pow(-0.0, ...) — знак результата. */
  ck_assert(s21_pow(-0.0, -3.0) == -INFINITY);
  ck_assert(s21_pow(-0.0, -2.0) == INFINITY);
  ck_assert(s21_pow(-0.0, 3.0) == 0.0);
  ck_assert(s21_pow(-0.0, 2.0) == 0.0);
}
END_TEST

START_TEST(test_pow_neg_large_exp) {
  ck_assert_ldouble_eq_tol(s21_pow(-2.0, 31.0), -2147483648.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(-2.0, 32.0), 4294967296.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(-2.0, 60.0), 1152921504606846976.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(-2.0, 61.0), -2305843009213693952.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(-1.0, 1000.0), 1.0, EPS);
  ck_assert_ldouble_eq_tol(s21_pow(-1.0, 1001.0), -1.0, EPS);
}
END_TEST

/* ---------- sin ---------- */
START_TEST(test_sin_loop) {
  double x = _i * 0.5 - 10.0;
  ck_assert_ldouble_eq_tol(s21_sin(x), sin(x), EPS);
}
END_TEST

START_TEST(test_sin_edge) {
  ck_assert_ldouble_eq_tol(s21_sin(0.0), 0.0, EPS);
  ck_assert_ldouble_eq_tol(s21_sin(S21_PI), sin(S21_PI), EPS);
  ck_assert_ldouble_eq_tol(s21_sin(S21_PI_2), 1.0, EPS);
  ck_assert_ldouble_eq_tol(s21_sin(1000.0), sin(1000.0), EPS);
  ck_assert(s21_sin(INFINITY) != s21_sin(INFINITY));
}
END_TEST

/* ---------- cos ---------- */
START_TEST(test_cos_loop) {
  double x = _i * 0.5 - 10.0;
  ck_assert_ldouble_eq_tol(s21_cos(x), cos(x), EPS);
}
END_TEST

START_TEST(test_cos_edge) {
  ck_assert_ldouble_eq_tol(s21_cos(0.0), 1.0, EPS);
  ck_assert_ldouble_eq_tol(s21_cos(S21_PI), -1.0, EPS);
  ck_assert_ldouble_eq_tol(s21_cos(S21_PI_2), cos(S21_PI_2), EPS);
  ck_assert_ldouble_eq_tol(s21_cos(1000.0), cos(1000.0), EPS);
  ck_assert(s21_cos(INFINITY) != s21_cos(INFINITY));
}
END_TEST

/* ---------- tan ---------- */
START_TEST(test_tan_loop) {
  double x = _i * 0.3 - 3.0;
  ck_assert_ldouble_eq_tol(s21_tan(x), tan(x), EPS);
}
END_TEST

START_TEST(test_tan_edge) {
  ck_assert_ldouble_eq_tol(s21_tan(0.0), 0.0, EPS);
  ck_assert_ldouble_eq_tol(s21_tan(S21_PI / 4.0), 1.0, EPS);
  ck_assert(s21_tan(INFINITY) != s21_tan(INFINITY));
}
END_TEST

/* ---------- atan ---------- */
START_TEST(test_atan_loop) {
  double x = _i * 0.5 - 10.0;
  ck_assert_ldouble_eq_tol(s21_atan(x), atan(x), EPS);
}
END_TEST

START_TEST(test_atan_edge) {
  ck_assert_ldouble_eq_tol(s21_atan(0.0), 0.0, EPS);
  ck_assert_ldouble_eq_tol(s21_atan(1.0), S21_PI / 4.0, EPS);
  ck_assert_ldouble_eq_tol(s21_atan(-1.0), -S21_PI / 4.0, EPS);
  ck_assert_ldouble_eq_tol(s21_atan(1000.0), atan(1000.0), EPS);
  ck_assert_ldouble_eq_tol(s21_atan(INFINITY), S21_PI_2, EPS);
  ck_assert_ldouble_eq_tol(s21_atan(-INFINITY), -S21_PI_2, EPS);
}
END_TEST

/* ---------- asin ---------- */
START_TEST(test_asin_loop) {
  double x = _i * 0.1 - 1.0;
  ck_assert_ldouble_eq_tol(s21_asin(x), asin(x), EPS);
}
END_TEST

START_TEST(test_asin_edge) {
  ck_assert_ldouble_eq_tol(s21_asin(0.0), 0.0, EPS);
  ck_assert_ldouble_eq_tol(s21_asin(1.0), S21_PI_2, EPS);
  ck_assert_ldouble_eq_tol(s21_asin(-1.0), -S21_PI_2, EPS);
  ck_assert(s21_asin(2.0) != s21_asin(2.0));
  ck_assert(s21_asin(-2.0) != s21_asin(-2.0));
}
END_TEST

/* ---------- acos ---------- */
START_TEST(test_acos_loop) {
  double x = _i * 0.1 - 1.0;
  ck_assert_ldouble_eq_tol(s21_acos(x), acos(x), EPS);
}
END_TEST

START_TEST(test_acos_edge) {
  ck_assert_ldouble_eq_tol(s21_acos(1.0), 0.0, EPS);
  ck_assert_ldouble_eq_tol(s21_acos(-1.0), S21_PI, EPS);
  ck_assert_ldouble_eq_tol(s21_acos(0.0), S21_PI_2, EPS);
  ck_assert(s21_acos(2.0) != s21_acos(2.0));
  ck_assert(s21_acos(-2.0) != s21_acos(-2.0));
}
END_TEST

/* ---------- ULP helper ---------- */

static long double s21_ulp_diff(long double got, long double expected) {
  if (got == expected) return 0.0L;
  if (got != got || expected != expected) return got - expected;

  long double ulp = S21_NEXTAFTER(expected, INFINITY) - expected;
  if (ulp <= 0.0L || ulp != ulp) {
    ulp = expected - S21_NEXTAFTER(expected, -INFINITY);
  }
  if (ulp <= 0.0L || ulp != ulp) {
    return fabsl(got - expected);
  }
  return fabsl(got - expected) / ulp;
}

/* ---------- precision comparison ---------- */
START_TEST(test_precision_compare) {
  const double args[] = {0.5,  1.0,  1.5,   2.0, 3.0, 5.0,
                         10.0, 20.0, 100.0, 1e3, 1e6, 1e10};
  const int n = (int)(sizeof(args) / sizeof(args[0]));

#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 113
  const char *mode = "128-bit long double (quad precision)";
#elif defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  const char *mode = "64-bit long double (double)";
#else
  const char *mode = "80-bit long double (x86 extended)";
#endif

  printf("\n=== Absolute deviation: |libm - s21_*| ===\n");
  printf("Mode: %s\n", mode);

  printf("%-8s", "x");
  for (int i = 0; i < n; i++) printf("%-11.2e", args[i]);
  printf("\n");

  printf("%-8s", "exp");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(S21_ABS_DIFF(S21_LIBM_FN(exp, args[i]), s21_exp(args[i])));
  printf("\n");

  printf("%-8s", "log");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(S21_ABS_DIFF(S21_LIBM_FN(log, args[i]), s21_log(args[i])));
  printf("\n");

  printf("%-8s", "sin");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(S21_ABS_DIFF(S21_LIBM_FN(sin, args[i]), s21_sin(args[i])));
  printf("\n");

  printf("%-8s", "cos");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(S21_ABS_DIFF(S21_LIBM_FN(cos, args[i]), s21_cos(args[i])));
  printf("\n");

  printf("%-8s", "atan");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(S21_ABS_DIFF(S21_LIBM_FN(atan, args[i]), s21_atan(args[i])));
  printf("\n");

  printf("\n=== ULP deviation: |libm - s21_*| / ulp(libm) ===\n");
  printf("Mode: %s\n", mode);
  printf("0 = exact, 1 = one ULP, 2 = two ULP, ...\n");

  printf("%-8s", "x");
  for (int i = 0; i < n; i++) printf("%-11.2e", args[i]);
  printf("\n");

  printf("%-8s", "exp");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(s21_ulp_diff(s21_exp(args[i]), S21_LIBM_FN(exp, args[i])));
  printf("\n");

  printf("%-8s", "log");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(s21_ulp_diff(s21_log(args[i]), S21_LIBM_FN(log, args[i])));
  printf("\n");

  printf("%-8s", "sin");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(s21_ulp_diff(s21_sin(args[i]), S21_LIBM_FN(sin, args[i])));
  printf("\n");

  printf("%-8s", "cos");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(s21_ulp_diff(s21_cos(args[i]), S21_LIBM_FN(cos, args[i])));
  printf("\n");

  printf("%-8s", "atan");
  for (int i = 0; i < n; i++)
    S21_PRINT_LD(s21_ulp_diff(s21_atan(args[i]), S21_LIBM_FN(atan, args[i])));
  printf("\n");

  /* Формальные проверки — допуск EPS, адаптированный под платформу */
  ck_assert_ldouble_eq_tol(s21_exp(1.0), S21_LIBM_FN(exp, 1.0), EPS);
  ck_assert_ldouble_eq_tol(s21_log(2.0), S21_LIBM_FN(log, 2.0), EPS);
  ck_assert_ldouble_eq_tol(s21_sin(1.0), S21_LIBM_FN(sin, 1.0), EPS);
  ck_assert_ldouble_eq_tol(s21_cos(1.0), S21_LIBM_FN(cos, 1.0), EPS);
  ck_assert_ldouble_eq_tol(s21_atan(1.0), S21_LIBM_FN(atan, 1.0), EPS);
}
END_TEST

/* ---------- suite ---------- */
Suite *s21_math_suite(void) {
  Suite *s = suite_create("s21_math");
  TCase *tc = tcase_create("core");

  tcase_add_test(tc, test_abs_basic);
  tcase_add_loop_test(tc, test_fabs_loop, 0, 40);
  tcase_add_loop_test(tc, test_ceil_loop, 0, 40);
  tcase_add_test(tc, test_ceil_edge);
  tcase_add_loop_test(tc, test_floor_loop, 0, 40);
  tcase_add_test(tc, test_floor_edge);
  tcase_add_loop_test(tc, test_fmod_loop, 0, 25);
  tcase_add_test(tc, test_fmod_edge);
  tcase_add_test(tc, test_fmod_large);
  tcase_add_loop_test(tc, test_sqrt_loop, 0, 30);
  tcase_add_test(tc, test_sqrt_edge);
  tcase_add_loop_test(tc, test_exp_loop, 0, 30);
  tcase_add_test(tc, test_exp_edge);
  tcase_add_loop_test(tc, test_log_loop, 0, 30);
  tcase_add_test(tc, test_log_edge);
  tcase_add_loop_test(tc, test_pow_loop, 0, 20);
  tcase_add_test(tc, test_pow_edge);
  tcase_add_test(tc, test_pow_neg_large_exp);
  tcase_add_loop_test(tc, test_sin_loop, 0, 40);
  tcase_add_test(tc, test_sin_edge);
  tcase_add_loop_test(tc, test_cos_loop, 0, 40);
  tcase_add_test(tc, test_cos_edge);
  tcase_add_loop_test(tc, test_tan_loop, 0, 20);
  tcase_add_test(tc, test_tan_edge);
  tcase_add_loop_test(tc, test_atan_loop, 0, 40);
  tcase_add_test(tc, test_atan_edge);
  tcase_add_loop_test(tc, test_asin_loop, 0, 21);
  tcase_add_test(tc, test_asin_edge);
  tcase_add_loop_test(tc, test_acos_loop, 0, 21);
  tcase_add_test(tc, test_acos_edge);
  tcase_add_test(tc, test_precision_compare);

  suite_add_tcase(s, tc);
  return s;
}

int main(void) {
  int failed = 0;
  Suite *s = s21_math_suite();
  SRunner *sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return failed == 0 ? 0 : 1;
}
