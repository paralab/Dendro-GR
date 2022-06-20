/**
 * @file dendro.cpp
 * @author Milinda Fernando
 * @brief dendro global definitions. 
 * @version 0.1
 * @date 2019-10-16
 * 
 * @copyright School of Computing, University of Utah. Copyright (c) 2019
 * 
 */
#include "dendro.h"

unsigned int MAXDEAPTH_LEVEL_DIFF =1;

void __handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);

}