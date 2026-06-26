CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRCS = src/circuit_elements.c src/dc_analysis.c src/ac_analysis.c src/transient.c src/circuit_matrix.c src/circuit_utils.c
OBJS = $(SRCS:.c=.o)
LIB = libminidcac.a

TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(TEST_SRCS:.c=.out)

EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=.out)

.PHONY: all test examples clean lib

all: lib test examples

lib: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

%.out: %.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lminidcac -o $@ $(LDFLAGS)

test: $(TEST_BINS)
	@echo "Running tests..."
	@for t in $(TEST_BINS); do echo "--- $$t ---"; ./$$t; done

examples: $(EXAMPLE_BINS)
	@echo "Running examples..."
	@for e in $(EXAMPLE_BINS); do echo "--- $$e ---"; ./$$e; done

clean:
	rm -f src/*.o $(LIB) tests/*.out examples/*.out demos/*.out benches/*.out

bench: benches/bench_matrix.out
	./benches/bench_matrix.out

demos: demos/demo_bode.out
	./demos/demo_bode.out
