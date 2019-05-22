#include "logger.h"
#include <cstdarg>
#include <cstdio>

PAT_log_level PAT_LOG_LEVEL = DEFAULT_LOGLEVEL;

char log_prefix[4][10] = {
    "[ERROR] ",
    "[WARN-] ",
    "[INFO-] ",
    "[DEBUG] ",
};


void log(PAT_log_level level, char *format, ...)
{
    if (level > PAT_LOG_LEVEL)
        return;

    printf("%s", log_prefix[level]);

    va_list args;

    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}
