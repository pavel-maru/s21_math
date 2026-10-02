CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11

# Опционально: 128-битный long double (только x86, медленнее).
# Обычно используется через цели test-128 / gcov-128 / clean-128.
# Можно также вручную:
#   make test LDBL=-mlong-double-128
#   make gcov_report LDBL=-mlong-double-128
LDBL ?=
CFLAGS += $(LDBL)

TEST_LIBS ?= -lcheck -lm -lpthread -lsubunit

SRC = s21_math.c
OBJ = s21_math.o
LIB = s21_math.a

FORMAT_FILES = s21_math.c s21_math.h test.c
CLANG_FORMAT ?= clang-format

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)

test: $(LIB)
	$(CC) $(CFLAGS) test.c $(LIB) $(TEST_LIBS) -o test
	./test

gcov_report: clean
	$(CC) $(CFLAGS) --coverage -c $(SRC) -o $(OBJ)
	$(CC) $(CFLAGS) --coverage test.c $(OBJ) $(TEST_LIBS) -o test
	./test
	lcov -t "s21_math" -o s21_math.info -c -d . --no-external
	genhtml -o report s21_math.info

# --- 128-битный long double (x86 only) ---
test-128:
	$(MAKE) test LDBL=-mlong-double-128

gcov-128:
	$(MAKE) gcov_report LDBL=-mlong-double-128

clean-128:
	$(MAKE) clean

# --- Форматирование ---
clang-format:
	$(CLANG_FORMAT) -i $(FORMAT_FILES)

clang-format-check:
	$(CLANG_FORMAT) --dry-run --Werror $(FORMAT_FILES)

# --- Проверки памяти (опционально) ---
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
	rm -f $(OBJ) $(LIB) test *.gcda *.gcno *.info
	rm -rf report

.PHONY: all clean test gcov_report test-128 gcov-128 clean-128 \
        clang-format clang-format-check valgrind asan ubsan
