#include "amatrix.h"
#include "global.h"

int amatrix_load(char *file, int a[])
{
    FILE *fp;
    char line[N+2];
    int i, j, k, m[N][N];

    if ((fp = fopen(file, "r")) == NULL) return 0;

    for (j = 0; j < N; j++)
    {
        if (fgets(line, N+2, fp) == NULL) return 0;

        for (k = 0; k < N; k++)
        {
            if (line[k] == '\0') return 0;
            m[j][k] = (line[k] == '0') ? 0 : 1;
        }
    }

    i = 0;

    for (j = 1; j < N; j++)
        for (k = 0; k < j; k++)
            a[i++] = m[j][k];

    return 1;
}
