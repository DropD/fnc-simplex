# Lightweight & simple Makefile - build 130529
# 2011 (C) Donjan Rodic - for free use
#
# Requires all source files & headers to be in one directory (SOURCEDIR)
# Produces one binary (BINDIR/BINARY)
#
# usage: make [ all | OBJ | debug | assembler | clean | install | uninstall
#                   | info ]
#  where OBJ is the path to an object file
# example: make
#          make install
#
################################################################################
# Config

# The directories (except SOURCEDIR & INFODIR) will be generated if not existing
SOURCEDIR = ./
OBJDIR = ../bin/objs
BINDIR = ../bin
BINARY = main
DESTDIR = ./
INFODIR = doc

# Path to headers for which (when changed) all objects should be recompiled
# Required if there is no source file with the same basename
SHAREDHEADERS = simplex/*.hpp misc/*.hpp *.hpp

# Only one SOURCEEXT is supported, but as many HEADEREXT as you wish
SOURCEEXT = cpp
HEADEREXT = h hpp hxx

CC = g++
CFLAGS = -Wall -Wextra -pedantic-errors -ansi
CFLAGS += -O3 -march=native -m64
CFLAGS += -fno-tree-vectorize -fno-tree-slp-vectorize -mno-abm

ifeq ($(shell uname), Darwin)
	CC = g++-4.7
endif

#~ CC = clang++
#~ CFLAGS = -pedantic -O3

CFLAGS += -std=c++11
CFLAGS += -I/opt/boost-1.53/include/
ifeq (!$(shell uname), Darwin)
	CFLAGS += -rdynamic
endif

CFLAGS += -I"../../gurobi/include"
CFLAGS += -I"../../soplex/src"

DBGFLAGS = $(CFLAGS) -Wall -Wextra -g3 -ggdb -pg
ifeq (!$(shell uname), Darwin)
	DBGFLAGS +=  -rdynamic
endif
LIBS = -lglpk -lm
LIBS += -L"../../gurobi/lib/" -lgurobi_c++ -lgurobi55 -lpthread
LIBS += -L"../../soplex/lib/" -lsoplex -lz

SMART_MKDIR = mkdir -p
SHELL = /bin/sh


################################################################################
# Automation starts here
# When trying to improve, don't remove "superfluous" code... it's needed!

.PHONY: all debug clean install uninstall info

BIN = $(BINDIR)/$(BINARY)
DESTBIN = $(DESTDIR)/$(BINARY)
VPATH = $(SOURCEDIR)
SOURCES=$(wildcard $(SOURCEDIR)/*.$(SOURCEEXT))
OBJS=$(subst $(SOURCEDIR),$(OBJDIR),$(SOURCES:.$(SOURCEEXT)=.o))
HEADERS:=$(foreach EXT,$(HEADEREXT), $(wildcard $(SOURCEDIR)/*.$(EXT)) )
HEADER_WITH_OBJ:=$(foreach FILE,$(HEADERS), \
$(if $(wildcard $(addsuffix .$(SOURCEEXT),$(basename $(FILE)) )), $(FILE)) \
)

all: $(BIN)

$(BIN): $(SHAREDHEADERS) $(OBJS)
	@$(SMART_MKDIR) $(BINDIR)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJS) $(LIBS)
	@echo Shared headers are: $(SHAREDHEADERS)

$(OBJDIR)/%.o: %.$(SOURCEEXT) $(SHAREDHEADERS)
	@$(SMART_MKDIR) $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Workaround: adds trigger on headers modification (obj.o: obj.HEADEREXT)
$(foreach FILE,$(HEADER_WITH_OBJ),$(eval \
$(subst $(SOURCEDIR), $(OBJDIR), $(addsuffix .o,$(basename $(FILE)) )): $(FILE)\
))


debug: $(SOURCES) $(HEADERS) $(SHAREDHEADERS)
	$(CC) $(DBGFLAGS) -o $(BIN)-dbg $(SOURCES) $(LIBS)

assembler: $(SOURCES) $(HEADERS) $(SHAREDHEADERS)
	@$(SMART_MKDIR) $(BINDIR)
	$(CC) $(CFLAGS) -S -o $(BIN).s $(SOURCES) $(LIBS)

clean:
	rm -f $(BIN)
	rm -f $(BIN)-dbg
	rm -f $(BIN).s
	rm -f $(BINDIR)/gmon.out
	rm -f $(OBJS)
	rm -rf $(OBJDIR)

install: $(BIN)
	@$(SMART_MKDIR) $(DESTDIR)
	cp $(BIN) $(DESTBIN)

uninstall: $(DESTBIN)
	rm $(DESTBIN)

info:
	cd $(INFODIR) && $(MAKE)
