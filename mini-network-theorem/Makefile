# Makefile for mini-network-theorem
# Network Theorems in Circuit Analysis

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm

SRCDIR   = src
TESTDIR  = tests
EXAMDIR  = examples
BUILDDIR = build

SRCS     = src/mesh_nodal.c src/thevenin_norton.c src/two_port.c src/superposition.c src/power_transfer.c src/wye_delta.c src/mna_solver.c
OBJS     = $(patsubst src/%.c,$(BUILDDIR)/%.o,$(SRCS))

.PHONY: all test examples clean

all: $(BUILDDIR)/libnetworktheorem.a

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.o: src/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/libnetworktheorem.a: $(OBJS)
	ar rcs $@ $^

test: $(BUILDDIR)/libnetworktheorem.a
	$(CC) $(CFLAGS) tests/test_all.c -L$(BUILDDIR) -lnetworktheorem $(LDFLAGS) -o $(BUILDDIR)/test_all
	./$(BUILDDIR)/test_all

examples: $(BUILDDIR)/libnetworktheorem.a
	$(CC) $(CFLAGS) examples/example_thevenin.c -L$(BUILDDIR) -lnetworktheorem $(LDFLAGS) -o $(BUILDDIR)/example_thevenin
	$(CC) $(CFLAGS) examples/example_superposition.c -L$(BUILDDIR) -lnetworktheorem $(LDFLAGS) -o $(BUILDDIR)/example_superposition
	$(CC) $(CFLAGS) examples/example_bridge.c -L$(BUILDDIR) -lnetworktheorem $(LDFLAGS) -o $(BUILDDIR)/example_bridge

clean:
	rm -rf $(BUILDDIR)
