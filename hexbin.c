void bin2hex(int *bin, int n, char *hex)
{
    int d, i = (4 - n%4)%4;

    while (n > 0)
    {
        d = 0;

        for (; i < 4; i++)
        {
            d = (d << 1) | *(bin++);
            if (--n == 0) break;
        }

        *(hex++) = (d < 10) ? d+'0' : d+'a'-10;
        i = 0;
    }

    *hex = '\0';
}

void hex2bin(char *hex, int n, int *bin)
{
    int d, i;

    for (; n > 0; n--)
    {
        d = *(hex++);
        d = (d < '0'+10) ? d-'0' : d-'a'+10;

        for (i = 0; i < 4; i++)
        {
            *(bin++) = (d & 8) != 0;
            d <<= 1;
        }
    }
}
