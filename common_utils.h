#ifndef PD_ASSERT_H
#define PD_ASSERT_H

#include <mpi.h>
#include <cassert>

#ifdef DEBUG
    #define PDASSERT(cond) \
    do { \
        assert(cond); \
    } while(false)
#else
    #define PDASSERT(cond) static_cast<void>(0)
#endif

#define PI ((double) 3.1415926535897932384626433)

#define PDLN_ABS_TOLERANCE_LOW ((double) 1e-11)
#define PDLN_ABS_TOLERANCE     ((double) 1e-10) // if less than 1e-10, will three point in a line, if more than 1e-15, will not pass check
#define PDLN_ABS_TOLERANCE_HI  ((double) 1e-11) // normal grid less than 1e-11

#define PDLN_RELATIVE_TOLERANCE_LOW ((double) 1e-3)
#define PDLN_RELATIVE_TOLERANCE     ((double) 1e-4)
#define PDLN_RELATIVE_TOLERANCE_HI  ((double) 1e-5)

#define float_eq_low(a, b) (fabs(a - b) <= PDLN_ABS_TOLERANCE_LOW)
#define float_eq(a, b)     (fabs(a - b) <= PDLN_ABS_TOLERANCE)
#define float_eq_hi(a, b)  (fabs(a - b) <= PDLN_ABS_TOLERANCE_HI)

#define float_relative_eq_low(a, b) (fabs((double)a - (double)b)/(double)a <= PDLN_RELATIVE_TOLERANCE_LOW)
#define float_relative_eq(a, b)     (fabs((double)a - (double)b)/(double)a <= PDLN_RELATIVE_TOLERANCE)
#define float_relative_eq_hi(a, b)  (fabs((double)a - (double)b)/(double)a <= PDLN_RELATIVE_TOLERANCE_HI)

class PDLN_Timer
{
    public:
        double time;

        PDLN_Timer() : time(0){};
        double tick() {
            double t = MPI_Wtime();
            double result = t - time;
            time = t;
            return result;
        }
};

#endif
