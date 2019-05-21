#include "logger.h"
#include <cstdarg>
#include <cstdio>

char log_prefix[4][10] = {
    "[ERROR-] ",
    "[-WARN-] ",
    "[REPORT] ",
    "[DEBUG-] ",
};


void log(PAT_LOG_LEVEL level, char *format, ...)
{
    printf("%s", log_prefix[level]);

    va_list args;

    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}
