#include "benchmark_utils.h"

#include <cstdio>
#include <iostream>

#include <fcntl.h>
#include <unistd.h>

StdoutSilencer::StdoutSilencer()
{
   std::cout.flush();
   std::fflush(stdout);
   const int null = open("/dev/null", O_WRONLY | O_CLOEXEC); // NOLINT(cppcoreguidelines-pro-type-vararg)
   if (null < 0)
      return;
   m_saved = dup(STDOUT_FILENO);
   if (m_saved >= 0 && dup2(null, STDOUT_FILENO) < 0) {
      close(m_saved);
      m_saved = -1;
   }
   close(null);
}

StdoutSilencer::~StdoutSilencer()
{
   if (m_saved < 0)
      return;
   std::cout.flush();
   std::fflush(stdout);
   dup2(m_saved, STDOUT_FILENO);
   close(m_saved);
}
