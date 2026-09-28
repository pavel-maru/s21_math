CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11

# Опционально: 128-битный long double (только x86, медленнее).
# Включается так:
#   make LDBL=-mlong-double-128
#   make test LDBL=-mlong-double-128
#   make gcov_report LDBL=-mlong-double-128
# Осторожно: libcheck может быть собран под 80-битный long double —
# тогда тесты начнут падать на сравнениях с эталоном. Проверяй перед сдачей.
LDBL ?=
CFLAGS += $(LDBL)

TEST_LIBS ?= -lcheck -lm -lpthread -lsubunit

SRC = s21_math.c
OBJ = s21_math.o
LIB = s21_math.a

# Файлы, к которым применяется clang-format
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

# Форматирование файлов по Google Style (см. .clang-format в корне src/)
clang-format:
	$(CLANG_FORMAT) -i $(FORMAT_FILES)

# Проверка без изменения файлов — удобно для CI и перед коммитом.
# Выход с ошибкой, если что-то не отформатировано.
clang-format-check:
	$(CLANG_FORMAT) --dry-run --Werror $(FORMAT_FILES)

clean:
	rm -f $(OBJ) $(LIB) test *.gcda *.gcno *.info
	rm -rf report

.PHONY: all clean test gcov_report clang-format clang-format-check
