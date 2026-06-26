# Makefile for mini-power-factor
# Power Factor Analysis Library — AC Power Theory
#
# References:
#   IEEE Std 1459-2010, IEEE Std 519-2014
#   Steinmetz (1893), Fortescue (1918), Park (1929)
#   Erickson & Maksimovic "Power Electronics" (2001)

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXM_DIR  = examples
BLD_DIR  = build

SRCS = $(SRC_DIR)/power_factor.c \
       $(SRC_DIR)/complex_power.c \
       $(SRC_DIR)/pf_correction.c \
       $(SRC_DIR)/harmonic_power.c \
       $(SRC_DIR)/power_quality.c \
       $(SRC_DIR)/phasor_operations.c

OBJS = $(BLD_DIR)/power_factor.o \
       $(BLD_DIR)/complex_power.o \
       $(BLD_DIR)/pf_correction.o \
       $(BLD_DIR)/harmonic_power.o \
       $(BLD_DIR)/power_quality.o \
       $(BLD_DIR)/phasor_operations.o

HDRS = $(INC_DIR)/power_factor.h \
       $(INC_DIR)/complex_power.h \
       $(INC_DIR)/pf_correction.h \
       $(INC_DIR)/harmonic_power.h \
       $(INC_DIR)/power_quality.h \
       $(INC_DIR)/phasor_operations.h

TEST_BIN  = $(BLD_DIR)/test_power_factor
EX1_BIN   = $(BLD_DIR)/example_industrial_pfc
EX2_BIN   = $(BLD_DIR)/example_harmonic_analysis
EX3_BIN   = $(BLD_DIR)/example_three_phase_pq

.PHONY: all clean test examples

all: $(BLD_DIR)/libpower_factor.a

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/libpower_factor.a: $(OBJS)
	ar rcs $@ $^

$(BLD_DIR)/power_factor.o: $(SRC_DIR)/power_factor.c $(HDRS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/complex_power.o: $(SRC_DIR)/complex_power.c $(INC_DIR)/complex_power.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/pf_correction.o: $(SRC_DIR)/pf_correction.c $(INC_DIR)/pf_correction.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/harmonic_power.o: $(SRC_DIR)/harmonic_power.c $(INC_DIR)/harmonic_power.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/power_quality.o: $(SRC_DIR)/power_quality.c $(INC_DIR)/power_quality.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/phasor_operations.o: $(SRC_DIR)/phasor_operations.c $(INC_DIR)/phasor_operations.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TST_DIR)/test_power_factor.c $(BLD_DIR)/libpower_factor.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lpower_factor $(LDFLAGS) -o $@

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(EX1_BIN): $(EXM_DIR)/example_industrial_pfc.c $(BLD_DIR)/libpower_factor.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lpower_factor $(LDFLAGS) -o $@

$(EX2_BIN): $(EXM_DIR)/example_harmonic_analysis.c $(BLD_DIR)/libpower_factor.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lpower_factor $(LDFLAGS) -o $@

$(EX3_BIN): $(EXM_DIR)/example_three_phase_pq.c $(BLD_DIR)/libpower_factor.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lpower_factor $(LDFLAGS) -o $@

clean:
	rm -rf $(BLD_DIR)
