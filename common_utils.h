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

#define PDLN_FLOAT_EQ_ERROR_LOW ((double) 1e-8)
#define PDLN_FLOAT_EQ_ERROR ((double) 1e-10) // if less than 1e-10, will three point in a line, if more than 1e-15, will not pass check
#define PDLN_FLOAT_EQ_ERROR_HI ((double) 1e-11) // normal grid less than 1e-11

#define float_eq(a, b)     (fabs(a - b) <= PDLN_FLOAT_EQ_ERROR)
#define float_eq_low(a, b) (fabs(a - b) <= PDLN_FLOAT_EQ_ERROR_LOW)
#define float_eq_hi(a, b)  (fabs(a - b) <= PDLN_FLOAT_EQ_ERROR_HI)

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
