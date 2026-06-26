# Makefile for mini-frequency-response
# Mini Electronic Info — Circuit Analysis: Frequency Response Module
#
# Targets:
#   make          — Build library (libfreqresp.a)
#   make test     — Build and run test suite
#   make examples — Build all examples
#   make clean    — Remove build artifacts

CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -O2 -g
CPPFLAGS = -Iinclude -D_GNU_SOURCE
LDFLAGS  = -lm

# Source files
SRC = src/frequency_response.c \
      src/transfer_function.c \
      src/bode_plot.c \
      src/filter_design.c \
      src/resonance.c \
      src/stability.c \
      src/network_function.c

# Object files
OBJ = $(SRC:.c=.o)

# Library
LIB = libfreqresp.a

# Test
TEST_SRC = tests/test_frequency_response.c
TEST_BIN = test_frequency_response

# Examples
EXAMPLES_SRC = examples/example_butterworth_filter.c \
               examples/example_rlc_resonance.c \
               examples/example_stability_analysis.c

EXAMPLES_BIN = $(EXAMPLES_SRC:.c=)

# Headers
HEADERS = include/frequency_response.h \
          include/transfer_function.h \
          include/bode_plot.h \
          include/filter_design.h \
          include/resonance.h \
          include/stability.h \
          include/network_function.h

.PHONY: all test examples clean

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

%.o: %.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Test target
test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_SRC) $(OBJ)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJ) $(LDFLAGS) -o $@

# Examples
examples: $(EXAMPLES_BIN)

examples/example_butterworth_filter: examples/example_butterworth_filter.c $(OBJ)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJ) $(LDFLAGS) -o $@

examples/example_rlc_resonance: examples/example_rlc_resonance.c $(OBJ)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJ) $(LDFLAGS) -o $@

examples/example_stability_analysis: examples/example_stability_analysis.c $(OBJ)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJ) $(LDFLAGS) -o $@

clean:
	rm -f src/*.o
	rm -f $(LIB)
	rm -f $(TEST_BIN)
	rm -f $(EXAMPLES_BIN)
	rm -f test_frequency_response.exe
	rm -f examples/*.exe
