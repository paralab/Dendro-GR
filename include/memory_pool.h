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

#ifndef MEMORY_POOL_H
#define MEMORY_POOL_H


#include <vector>

struct mem_blk {
public:
  unsigned int  size;
  bool          free;
  double*       data;
  
  mem_blk(unsigned int sz) {
    size = sz;
    free = true;
    data = NULL;
  }
  
  ~mem_blk() {
    if (data != NULL) { 
      delete [] data;
      data = NULL;
    }
  }
}; 

class memory_pool
{
private:
  unsigned long                       m_ulCount;
  
  unsigned int                        m_uiAlignment;
  unsigned int                        m_uiPadding;
  
  // for the entries ...
  // need, address, size, is_available. 
  std::vector<mem_blk>                m_vecBlocks; 
  
  
public:
  memory_pool(unsigned int pad = 3, unsigned int align = 32);
  
  ~memory_pool();
  
  // basic ops. block based only. of size (sz +2*pad)^3 
  double* allocate(unsigned int sz);
  void    free(double* buf);
  
  // how much memory is currently used by the pool
  unsigned long used();
  
  // delete all allocated blocks.
  void purge();
};

#endif // MEMORY_POOL_H
