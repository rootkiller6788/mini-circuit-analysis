# Makefile for mini-circuit-topology
# Circuit Topology Analysis Library

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EXM_DIR  = examples
BLD_DIR  = build

SRCS = $(SRC_DIR)/circuit_graph.c $(SRC_DIR)/circuit_mna.c $(SRC_DIR)/circuit_analysis.c $(SRC_DIR)/circuit_tellegen.c
OBJS = $(BLD_DIR)/circuit_graph.o $(BLD_DIR)/circuit_mna.o $(BLD_DIR)/circuit_analysis.o $(BLD_DIR)/circuit_tellegen.o
TEST_BIN = $(BLD_DIR)/test_topology
HDRS = $(INC_DIR)/circuit_topology.h $(INC_DIR)/circuit_graph.h $(INC_DIR)/circuit_mna.h $(INC_DIR)/circuit_analysis.h

.PHONY: all clean test examples

all: $(BLD_DIR)/libcircuit_topology.a

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/libcircuit_topology.a: $(OBJS)
	ar rcs $@ $^

$(BLD_DIR)/circuit_graph.o: $(SRC_DIR)/circuit_graph.c $(HDRS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/circuit_mna.o: $(SRC_DIR)/circuit_mna.c $(HDRS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/circuit_analysis.o: $(SRC_DIR)/circuit_analysis.c $(HDRS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(BLD_DIR)/circuit_tellegen.o: $(SRC_DIR)/circuit_tellegen.c $(HDRS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TST_DIR)/test_topology.c $(BLD_DIR)/libcircuit_topology.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lcircuit_topology $(LDFLAGS) -o $@

examples: $(BLD_DIR)/example_resistor_network $(BLD_DIR)/example_wheatstone_bridge $(BLD_DIR)/example_rlc_circuit

$(BLD_DIR)/example_resistor_network: $(EXM_DIR)/example_resistor_network.c $(BLD_DIR)/libcircuit_topology.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lcircuit_topology $(LDFLAGS) -o $@

$(BLD_DIR)/example_wheatstone_bridge: $(EXM_DIR)/example_wheatstone_bridge.c $(BLD_DIR)/libcircuit_topology.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lcircuit_topology $(LDFLAGS) -o $@

$(BLD_DIR)/example_rlc_circuit: $(EXM_DIR)/example_rlc_circuit.c $(BLD_DIR)/libcircuit_topology.a | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -L$(BLD_DIR) -lcircuit_topology $(LDFLAGS) -o $@

clean:
	rm -rf $(BLD_DIR)
