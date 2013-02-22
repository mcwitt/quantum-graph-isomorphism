#include "amatrix.h"

#define MAX_LINE_LEN 100

int amatrix_load(char *file, int num_vertices, int *a)
{
    FILE *fp;
    char line[MAX_LINE_LEN];
    int row, col;

    if ((fp = fopen(file, "r")) == NULL) return 0;

    for (row = 0; row < num_vertices; row++)
    {
        if (fgets(line, MAX_LINE_LEN, fp) == NULL) return 0;

        for (col = 0; col < num_vertices; col++)
        {
            if (line[col] == '\0') return 0;
            a[row * num_vertices + col] = (line[col] == '0') ? 0 : 1;
        }
    }

    return 1;
}
