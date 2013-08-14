#include "amatrix.h"

#define MAX_LINE_LEN 100

int amatrix_load(char *file, int n, int a[])
{
    FILE *fp;
    char line[MAX_LINE_LEN];
    int i, j, k, m[n][n];

    if ((fp = fopen(file, "r")) == NULL) return 0;

    for (j = 0; j < n; j++)
    {
        if (fgets(line, MAX_LINE_LEN, fp) == NULL) return 0;

        for (k = 0; k < n; k++)
        {
            if (line[k] == '\0') return 0;
            m[j][k] = (line[k] == '0') ? 0 : 1;
        }
    }

    i = 0;

    for (j = 1; j < n; j++)
        for (k = 0; k < j; k++)
            a[i++] = m[j][k];

    return 1;
}
