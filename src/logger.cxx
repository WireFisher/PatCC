#include "logger.h"
#include <cstdarg>
#include <cstdio>

PAT_LOG_LEVEL CURRENT_LOG_LEVEL = LOG_ERROR;

char log_prefix[4][10] = {
    "[ERROR] ",
    "[WARN-] ",
    "[INFO-] ",
    "[DEBUG] ",
};


void log(PAT_LOG_LEVEL level, char *format, ...)
{
    if (level > CURRENT_LOG_LEVEL)
        return;

    printf("%s", log_prefix[level]);

    va_list args;

    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}
