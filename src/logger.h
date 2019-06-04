/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef PAT_LOGGER
#define PAT_LOGGER

#include <cstdarg>

#ifndef DEFAULT_LOGLEVEL
#define DEFAULT_LOGLEVEL LOG_INFO
#endif

extern char log_prefix[4][10];

enum PAT_log_level {
    LOG_ERROR = 0,
    LOG_WARNING = 1,
    LOG_INFO = 2,
    LOG_DEBUG = 3,
    LOG_DEBUG_V = 4,
};

extern PAT_log_level CURRENT_LOG_LEVEL;
extern void log(PAT_log_level level, char *format, ...);

#endif
