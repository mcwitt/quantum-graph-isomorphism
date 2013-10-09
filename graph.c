#include "graph.h"
#include "hexbin.h"
#include <math.h>

void graph_to_hex(graph_t *g, char *hex)
{
    bin2hex(g->b, N*(N-1)/2, hex);
}

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

void graph_read_hex(graph_t *g, char *hex)
{
    hex2bin(hex, (int) ceil(N*(N-1)/8), g->b);
}
