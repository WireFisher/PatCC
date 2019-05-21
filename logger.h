#ifndef PAT_LOGGER
#define PAT_LOGGER

#include <cstdarg>

extern char log_prefix[4][10];

enum PAT_LOG_LEVEL {
    LOG_ERROR = 0,
    LOG_WARNING = 1,
    LOG_REPORT = 2,
    LOG_DEBUG = 3,
};

extern void log(PAT_LOG_LEVEL level, char *format, ...);

#endif
