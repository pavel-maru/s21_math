CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11
# Опционально: включает 128-битный long double (x86 only, может замедлить тесты)
# CFLAGS += -mlong-double-128

TEST_LIBS = -lcheck -lm -lpthread -lsubunit

SRC = s21_math.c
OBJ = s21_math.o
LIB = s21_math.a

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

clean:
	rm -f $(OBJ) $(LIB) test *.gcda *.gcno *.info
	rm -rf report
