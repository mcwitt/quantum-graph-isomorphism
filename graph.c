#include "graph.h"
#include "hexbin.h"
#include <math.h>
#include <string.h>

int graph_read_amatrix(graph_t *g, FILE *fp)
{
    char line[N+2];
    int i, j, k, m[N][N];

    for (j = 0; j < N; j++)
    {
        if (fgets(line, N+2, fp) == NULL) return 1;

        for (k = 0; k < N; k++)
        {
            if (line[k] == '\0') return 1;
            m[j][k] = (line[k] == '0') ? 0 : 1;
        }
    }

    do fgets(line, N+2, fp); while (line[0] == '\n');
    if (! feof(fp)) return 1;

    i = 0;

    for (j = 1; j < N; j++)
        for (k = 0; k < j; k++)
            g->b[i++] = m[j][k];

    return 0;
}

int graph_read_bits(graph_t *g, char *bits)
{
    int i;

    if (strlen(bits) != GRAPH_BITS_LEN) return 1;

    for (i = 0; i < GRAPH_BITS_LEN; i++)
        g->b[i] = (bits[i] == '0') ? 0 : 1;

    return 0;
}

int graph_read_hexs(graph_t *g, char *hexs)
{
    int hexs_len = ceil(GRAPH_BITS_LEN / 4.);
    if (strlen(hexs) != hexs_len) return 1;
    return hex2bin(hexs, hexs_len, g->b);
}

void graph_to_hexs(graph_t *g, char *hexs)
{
    bin2hex(g->b, N*(N-1)/2, hexs);
}
