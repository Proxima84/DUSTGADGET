#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file system.c
 *  \brief contains miscellaneous routines, e.g. elapsed time measurements
 */

/*! returns the maximum of two double
 */
double dmax(double x, double y)
{
    if (x > y)
        return x;
    else
        return y;
}

/*! returns the minimum of two double
 */
double dmin(double x, double y)
{
    if (x < y)
        return x;
    else
        return y;
}

/*! returns the maximum of two integers
 */
int imax(int x, int y)
{
    if (x > y)
        return x;
    else
        return y;
}

/*! returns the minimum of two integers
 */
int imin(int x, int y)
{
    if (x < y)
        return x;
    else
        return y;
}
