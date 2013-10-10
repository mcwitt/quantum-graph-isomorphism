#ifndef GRAPH_H
#define GRAPH_H

#include "global.h"
#include <stdio.h>

#define GRAPH_BITS_LEN N*(N-1)/2

/**
 * Representation of an undirected graph.
 */
typedef struct
{
    int b[GRAPH_BITS_LEN];  /* independent adjacency matrix entries
                            (A_21, A_31, A_32, A_41, etc.) */
} graph_t;

/**
 * Load adjacency matrix from a text file
 * @param[out]  g       Graph struct.
 * @param[in]   file    Input file.
 */
int graph_read_amatrix(graph_t *g, FILE *fp);
int graph_read_bits(graph_t *g, char *bits);
int graph_read_hexs(graph_t *g, char *hexs);
void graph_to_hexs(graph_t *g, char *hexs);

#endif
