#pragma once
#include <cassert>
double Q[1220];
int    indx;
double cc;
double c;    /* current CSWB */
double zc;   /* current SWB `borrow` */
double zx;   /* SWB seed1 */
double zy;   /* SWB seed2 */
size_t qlen; /* length of Q array */

inline void InitiateRandNumGenerator(unsigned long x)
{


    assert(sizeof(double) >= 8);
    cc = 1.0 / 9007199254740992.0;  // inverse of 2^53rd power
    int    i;
    size_t qlen = indx = sizeof Q / sizeof Q[0];
    for (i = 0; i < qlen; i++)
        Q[i] = 0;
    double c = 0.0, zc = 0.0, /* current CSWB and SWB `borrow` */
        zx = 5212886298506819.0 / 9007199254740992.0, /* SWB seed1 */
        zy = 2020898595989513.0 / 9007199254740992.0; /* SWB seed2 */
    int                              j;
    double                           s, t; /* Choose 32 bits for x, 32 for y */
    unsigned long /*x = 123456789,*/ y =
        362436069; /* default seeds * /
               /* Next, seed each Q[i], one bit at a time, */

    if (x == 0) {
        x = 123456789;
    }

    for (i = 0; i < qlen; i++) { /* using 9th bit from Cong+Xorshift */
        s = 0.0;
        t = 1.0;
        for (j = 0; j < 52; j++) {
            t = 0.5 * t; /* make t=.5/2^j */
            x = 69069 * x + 123;
            y ^= (y << 13);
            y ^= (y >> 17);
            y ^= (y << 5);
            if (((x + y) >> 23) & 1)
                s = s + t; /* change bit of s, maybe */
        }                  /* end j loop */
        Q[i] = s;
    } /* end i seed loop, Now generate 10^9 RandNumGenerator's: */
}

inline double RandNumGenerator()
{ /* Takes 14 nanosecs, Intel Q6600,2.40GHz */
    int    i, j;
    double t; /* t: first temp, then next CSWB value */
    /* First get zy as next lag-2 SWB */
    t = zx - zy - zc;
    zx = zy;
    if (t < 0) {
        zy = t + 1.0;
        zc = cc;
    } else {
        zy = t;
        zc = 0.0;
    }

    /* Then get t as the next lag-1220 CSWB value */
    if (indx < 1220)
        t = Q[indx++];
    else { /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
        for (i = 0; i < 1220; i++) {
            j = (i < 30) ? i + 1190 : i - 30;
            t = Q[j] - Q[i] + c; /* Get next CSWB element */
            if (t > 0) {
                t = t - cc;
                c = cc;
            } else {
                t = t - cc + 1.0;
                c = 0.0;
            }
            Q[i] = t;
        } /* end i loop */
        indx = 1;
        t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
    }             /* end else segment; return t-zy mod 1 */
    return ((t < zy) ? 1.0 + (t - zy) : t - zy);
} /* end RandNumGenerator() */