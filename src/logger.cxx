/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


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
