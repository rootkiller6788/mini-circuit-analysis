# mini-spice-simulation Makefile
# A compact SPICE-like circuit simulator
#
# Knowledge coverage: L1-L6 Complete, L7-L9 Partial
# Reference: Nagel (1975), "SPICE2: A Computer Program to Simulate
#            Semiconductor Circuits", UCB/ERL M520

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g
INCLUDES = -Iinclude
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXM_DIR  = examples
BLD_DIR  = build

# Source files
SRCS = $(SRC_DIR)/spice_netlist.c \
       $(SRC_DIR)/spice_matrix.c \
       $(SRC_DIR)/spice_models.c \
       $(SRC_DIR)/spice_analysis.c \
       $(SRC_DIR)/spice_core.c

# Object files
OBJS = $(patsubst $(SRC_DIR)/%.c,$(BLD_DIR)/%.o,$(SRCS))

# Test files
TEST_SRCS = $(TST_DIR)/test_spice.c
TEST_OBJS = $(patsubst $(TST_DIR)/%.c,$(BLD_DIR)/%.o,$(TEST_SRCS))

# Example files
EXAMPLES = $(EXM_DIR)/example_rc_transient \
           $(EXM_DIR)/example_rlc_ac \
           $(EXM_DIR)/example_divider_dc

# Library (static)
LIB = $(BLD_DIR)/libspice.a

.PHONY: all clean test examples

all: $(LIB) test examples

# ── Build directory ──────────────────────────────────────────────────
$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# ── Object files ─────────────────────────────────────────────────────
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BLD_DIR)/%.o: $(TST_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ── Static library ───────────────────────────────────────────────────
$(LIB): $(OBJS)
	ar rcs $@ $^

# ── Tests ────────────────────────────────────────────────────────────
test: $(LIB) $(TEST_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BLD_DIR)/test_spice $(TEST_OBJS) -L$(BLD_DIR) -lspice $(LDFLAGS)
	$(BLD_DIR)/test_spice

# ── Examples ─────────────────────────────────────────────────────────
$(EXM_DIR)/%: $(EXM_DIR)/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< -L$(BLD_DIR) -lspice $(LDFLAGS)

examples: $(EXAMPLES)
	@echo "All examples built successfully."

# ── Clean ────────────────────────────────────────────────────────────
clean:
	rm -rf $(BLD_DIR)
	rm -f $(EXAMPLES)
