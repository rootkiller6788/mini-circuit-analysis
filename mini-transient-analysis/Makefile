CC = gcc
CFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -std=gnu11 -O2 -Iinclude
LDFLAGS = -lm

SRCS = src/transient_defs.c src/first_order.c src/second_order.c src/higher_order.c src/numerical_transient.c src/switching_transient.c
OBJS = $(SRCS:.c=.o)
LIB = libminitransient.a

TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(TEST_SRCS:.c=.out)

EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=.out)

.PHONY: all test examples clean lib bench demos

all: lib test examples

lib: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

%.out: %.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lminitransient -o $@ $(LDFLAGS)

test: $(TEST_BINS)
	@echo "Running tests..."
	@for t in $(TEST_BINS); do echo "--- $$t ---"; ./$$t; done

examples: $(EXAMPLE_BINS)
	@echo "Running examples..."
	@for e in $(EXAMPLE_BINS); do echo "--- $$e ---"; ./$$e; done

clean:
	rm -f src/*.o $(LIB) tests/*.out examples/*.out demos/*.out benches/*.out

bench: benches/bench_transient.out
	./benches/bench_transient.out

demos: demos/demo_rc_discharge.out
	./demos/demo_rc_discharge.out
