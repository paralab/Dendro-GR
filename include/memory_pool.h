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

#pragma once
#include <vector>

namespace mem
{
  /**
   * @brief simple structure to keep the mem blocks
   * 
   * @tparam T 
   */
  template<typename T>
  struct mem_blk {
    
    public:
      /**@brief : size of the block */
      unsigned int  size;
      /**@brief: true if currently in use false otherwise*/
      bool          free;
      /**@brief: pointer to the current block. */
      T*       data;
  
    mem_blk(unsigned int sz) {
      size = sz;
      free = true;
      data = NULL;
    }
  
    ~mem_blk() {}
  };

  template<typename T>
  class memory_pool
  {
    private:
      unsigned long                       m_ulCount;
      unsigned int                        m_uiAlignment;
      unsigned int                        m_uiPadding;
      std::vector<mem_blk<T>>             m_vecBlocks; 
    
    public:
      memory_pool(unsigned int pad = 3, unsigned int align = 32)
      {
        m_uiPadding = pad;
        m_uiAlignment = align;

      }

      ~memory_pool()
      {
        this->purge();        
      }
    
    /**
     * @brief allocate memory
     * 
     * @param sz 
     * @return T* 
     */
    T* allocate(unsigned int sz)
    {

      #if 1
        // if want allocatations without the memeory pool. 
        return new T[sz];
      #endif

      for(unsigned int i=0; i < m_vecBlocks.size(); i++) {
        if (m_vecBlocks[i].free && m_vecBlocks[i].size == sz) {
          m_vecBlocks[i].free = false;
          return m_vecBlocks[i].data;
        }
      }
      
      // 2. did not find a free block, create a new one.
      mem_blk<T> blk(sz);
      T *buf; 
      unsigned long len = sz;//(sz+m_uiPadding)*(sz+m_uiPadding)*(sz+m_uiPadding);
  
      // todo @hari @milinda allocate properly to be aligned in 3d
      const int st=posix_memalign((void **)&buf, m_uiAlignment, (len)*sizeof(T)); 
      blk.data = buf;
      blk.free = false;
      m_vecBlocks.push_back(blk);
      m_ulCount += len*sizeof(T);
      return buf;
    }

    void free(T* buf)
    {

      #if 1 
        delete [] buf;
        return;
      #endif

      for(unsigned int i=0; i < m_vecBlocks.size(); i++) {
        if (m_vecBlocks[i].data == buf) {
          m_vecBlocks[i].free = true;
          unsigned long len = m_vecBlocks[i].size;//(x.size+m_uiPadding)*(x.size+m_uiPadding)*(x.size+m_uiPadding);
          m_ulCount -= len*sizeof(T);
          return;
        }
      }
     
      std::cout << "memory_pool error: trying to free block not allocated by mempool."<<std::endl;
    }
    
    /**@brief how much memory is currently used by the pool*/
    unsigned long used(){ return m_ulCount;}
    
    /**
     * @brief delete all allocated blocks.
     */
    void purge()
    {
      for(unsigned int i=0; i < m_vecBlocks.size(); i++) {
          if (m_vecBlocks[i].data != NULL) {
            delete [] m_vecBlocks[i].data;
            m_vecBlocks[i].data = NULL;
          }
        }
      m_vecBlocks.clear();
    }
  };

} // end of namespace mem. 