#ifndef PD_ASSERT_H
#define PD_ASSERT_H

#include <cassert>

#ifdef DEBUG
    #define PDASSERT(cond) \
    do { \
        assert(cond); \
    } while(false)
#else
    #define PDASSERT(cond) static_cast<void>(0)
#endif

#endif
