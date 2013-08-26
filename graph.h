#include "global.h"
#include <stdio.h>

/**
 * Representation of an undirected graph.
 */
typedef struct
{
    int b[N*(N-1)/2];   /* independent adjacency matrix entries
                           (A_21, A_31, A_32, A_41, etc.) */
} graph_t;

/**
 * Load adjacency matrix from a text file
 * @param[out]  g       Graph struct.
 * @param[in]   file    Input file.
 */
int graph_read_amatrix(graph_t *g, char *file);

