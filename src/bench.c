#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "s21_math.h"

#define ITERS 1000000
#define NARGS 1000

static double pos_args[NARGS];  /* 0.1 .. 1.1 */
static double any_args[NARGS];  /* -5.0 .. 4.0 */
static double unit_args[NARGS]; /* -0.9 .. 0.9 */

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static volatile long double sink = 0.0L;

/* При -mlong-double-64 на x86-64 glibc собрана под 80-битный ABI.
   Прямой вызов expl/logl/sinl/... читает из регистров 16 байт и
   получает мусор. Используем double-версии и расширяем результат. */
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  #define S21_LIBM1(name, x)     ((long double)name((double)(x)))
  #define S21_LIBM2(name, x, y)  ((long double)name((double)(x), (double)(y)))
#else
  #define S21_LIBM1(name, x)     name##l(x)
  #define S21_LIBM2(name, x, y)  name##l((x), (y))
#endif

/* Обёртки: bench_unary/bench_binary принимают указатель на функцию,
   а макрос развернуть в указатель нельзя. */
static long double libm_exp (long double x) { return S21_LIBM1(exp,  x); }
static long double libm_log (long double x) { return S21_LIBM1(log,  x); }
static long double libm_sqrt(long double x) { return S21_LIBM1(sqrt, x); }
static long double libm_sin (long double x) { return S21_LIBM1(sin,  x); }
static long double libm_cos (long double x) { return S21_LIBM1(cos,  x); }
static long double libm_tan (long double x) { return S21_LIBM1(tan,  x); }
static long double libm_atan(long double x) { return S21_LIBM1(atan, x); }
static long double libm_asin(long double x) { return S21_LIBM1(asin, x); }
static long double libm_acos(long double x) { return S21_LIBM1(acos, x); }
static long double libm_pow (long double x, long double y) {
  return S21_LIBM2(pow, x, y);
}
static long double libm_fmod(long double x, long double y) {
  return S21_LIBM2(fmod, x, y);
}

typedef long double (*s21_unary_fn)(double);
typedef long double (*libm_unary_fn)(long double);

static void bench_unary(const char *name, s21_unary_fn s21_fn,
                        libm_unary_fn libm_fn, const double *args) {
  for (int i = 0; i < NARGS; i++) {
    sink += s21_fn(args[i]);
    sink += libm_fn((long double)args[i]);
  }

  double t0 = now_sec();
  for (int i = 0; i < ITERS; i++) {
    sink += s21_fn(args[i % NARGS]);
  }
  double t1 = now_sec();

  double t2 = now_sec();
  for (int i = 0; i < ITERS; i++) {
    sink += libm_fn((long double)args[i % NARGS]);
  }
  double t3 = now_sec();

  double s21_ns = (t1 - t0) * 1e9 / (double)ITERS;
  double libm_ns = (t3 - t2) * 1e9 / (double)ITERS;
  double ratio = (libm_ns > 0.0) ? s21_ns / libm_ns : 0.0;

  printf("%-8s %14.2f %14.2f %10.2fx\n", name, s21_ns, libm_ns, ratio);
}

typedef long double (*s21_binary_fn)(double, double);
typedef long double (*libm_binary_fn)(long double, long double);

static void bench_binary(const char *name, s21_binary_fn s21_fn,
                         libm_binary_fn libm_fn, const double *args,
                         double second) {
  for (int i = 0; i < NARGS; i++) {
    sink += s21_fn(args[i], second);
    sink += libm_fn((long double)args[i], (long double)second);
  }

  double t0 = now_sec();
  for (int i = 0; i < ITERS; i++) {
    sink += s21_fn(args[i % NARGS], second);
  }
  double t1 = now_sec();

  double t2 = now_sec();
  for (int i = 0; i < ITERS; i++) {
    sink += libm_fn((long double)args[i % NARGS], (long double)second);
  }
  double t3 = now_sec();

  double s21_ns = (t1 - t0) * 1e9 / (double)ITERS;
  double libm_ns = (t3 - t2) * 1e9 / (double)ITERS;
  double ratio = (libm_ns > 0.0) ? s21_ns / libm_ns : 0.0;

  printf("%-8s %14.2f %14.2f %10.2fx\n", name, s21_ns, libm_ns, ratio);
}

int main(void) {
  for (int i = 0; i < NARGS; i++) {
    pos_args[i] = 0.1 + (double)i * 0.001;
    any_args[i] = -5.0 + (double)i * 0.009;
    unit_args[i] = -0.9 + (double)i * 0.0018;
  }

  printf("Benchmark: %d calls per function\n", ITERS);
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 113
  printf("Mode: 128-bit long double (quad precision)\n\n");
#elif defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 53
  printf("Mode: 64-bit long double (double)\n\n");
#else
  printf("Mode: 80-bit long double (x86 extended)\n\n");
#endif

  printf("%-8s %14s %14s %11s\n", "fn", "s21 (ns/call)", "libm (ns/call)",
         "ratio");
  printf("--------------------------------------------------------------\n");

  /* Одноместные — положительные аргументы */
  bench_unary("exp",  s21_exp,  libm_exp,  pos_args);
  bench_unary("log",  s21_log,  libm_log,  pos_args);
  bench_unary("sqrt", s21_sqrt, libm_sqrt, pos_args);

  /* Одноместные — произвольные аргументы */
  bench_unary("sin",  s21_sin,  libm_sin,  any_args);
  bench_unary("cos",  s21_cos,  libm_cos,  any_args);
  bench_unary("tan",  s21_tan,  libm_tan,  any_args);
  bench_unary("atan", s21_atan, libm_atan, any_args);

  /* asin/acos — аргументы в области определения [-1, 1] */
  bench_unary("asin", s21_asin, libm_asin, unit_args);
  bench_unary("acos", s21_acos, libm_acos, unit_args);

  /* Бинарные */
  bench_binary("pow",  s21_pow,  libm_pow,  pos_args, 1.5);
  bench_binary("fmod", s21_fmod, libm_fmod, pos_args, 0.3);

  printf("--------------------------------------------------------------\n");
  printf("(ratio > 1 — s21_* медленнее libm, ratio < 1 — быстрее)\n");

  if (sink == 42.0L) printf("unreachable: %Lf\n", sink);

  return 0;
}
