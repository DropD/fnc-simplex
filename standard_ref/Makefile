.PHONY: all what

CC = g++
CFLAGS = -Wall -Wextra -Werror -pedantic-errors -ansi
CFLAGS += -Ofast -march=native -m64

CC = clang++
CFLAGS = -pedantic -O3


CFLAGS += -std=c++11
CFLAGS += -I/opt/boost-1.53/include/
CFLAGS += -rdynamic

BINDIR = ./

all:
	@echo "type 'make what' for targets"

%:	%.c*  *.hpp
	`which time` -f "\treal: %e \tuser: %U \tsys: %S  \twall clock: %e" \
	$(CC) $(CFLAGS) $< -o $(BINDIR)/$@


what:
	@echo "AVAILABLE TARGETS:"
	@ls | grep -E "(\.c|\.cpp)" | sed 's/\..*//'
