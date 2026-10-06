CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11

# Опционально: разрядность long double (только x86-64, медленнее).
#   LDBL=-mlong-double-128  — 128-бит quad (эмуляция ARM64 Linux)
#   LDBL=-mlong-double-64   — 64-бит double (эмуляция Apple Silicon)
LDBL ?=
CFLAGS += $(LDBL)

TEST_LIBS ?= -lcheck -lm -lpthread -lsubunit

SRC = s21_math.c
OBJ = s21_math.o
LIB = s21_math.a

FORMAT_FILES = s21_math.c s21_math.h test.c bench.c
CLANG_FORMAT ?= clang-format

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)

# test всегда пересобирает библиотеку с нуля. Это защищает от
# ситуации, когда s21_math.a остался от прогона в другой
# разрядности long double: несовпадение ABI даёт мусор в тестах.
test:
	$(MAKE) clean
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)
	ar rcs $(LIB) $(OBJ)
	$(CC) $(CFLAGS) test.c $(LIB) $(TEST_LIBS) -o test
	./test

gcov_report: clean
	$(CC) $(CFLAGS) --coverage -c $(SRC) -o $(OBJ)
	$(CC) $(CFLAGS) --coverage test.c $(OBJ) $(TEST_LIBS) -o test
	./test
	lcov -t "s21_math" -o s21_math.info -c -d . --no-external
	genhtml -o report s21_math.info

# --- 128-битный long double (только x86-64 Linux, gcc) ---
test-128:
	$(MAKE) test LDBL=-mlong-double-128

gcov-128:
	$(MAKE) gcov_report LDBL=-mlong-double-128

clean-128:
	$(MAKE) clean

# --- 64-битный long double (симуляция Apple Silicon) ---
test-64:
	$(MAKE) test LDBL=-mlong-double-64

gcov-64:
	$(MAKE) gcov_report LDBL=-mlong-double-64

clean-64:
	$(MAKE) clean

# --- Бенчмарк ---
bench: clean
	$(CC) $(CFLAGS) -O2 bench.c $(SRC) -lm -o bench
	./bench
	rm -f bench

# --- Форматирование ---
clang-format:
	$(CLANG_FORMAT) -i $(FORMAT_FILES)

clang-format-check:
	$(CLANG_FORMAT) --dry-run --Werror $(FORMAT_FILES)

# --- Проверки памяти ---
valgrind: test
	valgrind --tool=memcheck --leak-check=yes --show-leak-kinds=definite ./test

asan: clean
	$(CC) $(CFLAGS) -fsanitize=address -g -c $(SRC) -o $(OBJ)
	$(CC) $(CFLAGS) -fsanitize=address -g test.c $(OBJ) $(TEST_LIBS) -o test
	./test

ubsan: clean
	$(CC) $(CFLAGS) -fsanitize=undefined -g -c $(SRC) -o $(OBJ)
	$(CC) $(CFLAGS) -fsanitize=undefined -g test.c $(OBJ) $(TEST_LIBS) -o test
	./test

clean:
	rm -f $(OBJ) $(LIB) test bench *.gcda *.gcno *.info
	rm -rf report

.PHONY: all clean test gcov_report test-128 gcov-128 clean-128 \
        test-64 gcov-64 clean-64 \
        clang-format clang-format-check valgrind asan ubsan bench
