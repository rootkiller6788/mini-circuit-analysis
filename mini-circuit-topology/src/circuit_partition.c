/** @file circuit_partition.c
 * @brief Circuit partitioning for hierarchical analysis
 * Splits large circuits into subcircuits for parallel/distributed solving.
 * Knowledge: L5 (Algorithms): Kernighan-Lin partitioning
 *             L8 (Advanced): Hierarchical circuit simulation
 */
#include "circuit_topology.h"
#include "circuit_graph.h"
#include <string.h>
#include <stdlib.h>

int ct_partition_by_degree(const ct_circuit_t *circuit, int32_t *partition,
                           int32_t num_parts, int32_t *part_sizes)
{
    if (!circuit || !partition || !part_sizes || num_parts < 2) return -1;
    int32_t n = circuit->num_nodes;
    for (int32_t i = 0; i < num_parts; i++) part_sizes[i] = 0;
    partition[0] = 0; part_sizes[0]++;
    for (int32_t i = 1; i < n; i++) {
        int deg = ct_node_degree(circuit, i);
        int32_t best_part = 0;
        int32_t min_size = part_sizes[0];
        for (int32_t p = 0; p < num_parts; p++) {
            if (part_sizes[p] < min_size) { min_size = part_sizes[p]; best_part = p; }
        }
        partition[i] = best_part;
        part_sizes[best_part]++;
    }
    return 0;
}

int ct_count_cut_edges(const ct_circuit_t *circuit, const int32_t *partition,
                       int32_t num_parts, int32_t *cut_count)
{
    if (!circuit || !partition || !cut_count) return -1;
    *cut_count = 0;
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        int32_t nf = circuit->branches[b].node_from;
        int32_t nt = circuit->branches[b].node_to;
        if (partition[nf] != partition[nt]) (*cut_count)++;
    }
    return 0;
}

double ct_partition_balance(const int32_t *part_sizes, int32_t num_parts,
                           int32_t total_nodes)
{
    if (!part_sizes || num_parts < 1 || total_nodes <= 0) return -1.0;
    double avg = (double)total_nodes / num_parts;
    double max_dev = 0.0;
    for (int32_t i = 0; i < num_parts; i++) {
        double dev = (part_sizes[i] - avg) / avg;
        if (dev < 0) dev = -dev;
        if (dev > max_dev) max_dev = dev;
    }
    return 1.0 - max_dev;
}