/*
 * Copyright (c) 2017 Hari Sundar <hari@cs.utah.edu>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include "memory_pool.h"

#include <cstdlib>

memory_pool::memory_pool(unsigned int pad, unsigned int align) {
  m_uiPadding = pad;
  m_uiAlignment = align;
}

memory_pool::~memory_pool() {
  
}

// basic ops. block based only. of size (sz +2*pad)^3 
double* 
memory_pool::allocate(unsigned int sz) {
  // 1. find if a free block of size exists.
  for (auto x: m_vecBlocks) {
    if (x.size == sz && x.free) {
      x.free = false;
      return x.data;
    }
  }
  // 2. did not find a free block, create a new one.
  mem_blk blk(sz);
  double *buf; 
  unsigned long len = (sz+m_uiPadding)*(sz+m_uiPadding)*(sz+m_uiPadding);
  
  // todo @hari allocate properly to be aligned in 3d
  posix_memalign((void **)&buf, m_uiAlignment, (len)*sizeof(double)); 
  
  blk.data = buf;
  blk.free = false;
  
  m_vecBlocks.push_back(blk);
  
  m_ulCount += len*sizeof(double);
  
  return buf;
}

void    
memory_pool::free(double* buf) {
  for (auto x: m_vecBlocks) {
    if (x.data == buf) {
      x.free = true;
      unsigned long len = (x.size+m_uiPadding)*(x.size+m_uiPadding)*(x.size+m_uiPadding);
      m_ulCount -= len*sizeof(double);
      return;
    }
  }
  std::cout << "memory_pool error: trying to free block not allocated by mempool."
}

// how much memory is currently used by the pool
unsigned long 
memory_pool::used() {
  return m_ulCount;
}

// delete all allocated blocks.
void 
memory_pool::purge() {
  for (auto x: m_vecBlocks) {
    if (x.data != NULL) {
      delete [] x.data;
      x.data = NULL;
    }
  }
  m_vecBlocks.clear();
}
