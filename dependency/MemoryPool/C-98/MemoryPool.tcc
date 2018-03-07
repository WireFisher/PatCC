/*-
 * Copyright (c) 2013 Cosku Acay, http://www.coskuacay.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MEMORY_BLOCK_TCC
#define MEMORY_BLOCK_TCC

#include <cstdio>

//template <typename T, size_t BlockSize>
template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::size_type
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::padPointer(data_pointer_ p, size_type align)
const throw()
{
  size_t result = reinterpret_cast<size_t>(p);
  return ((align - result) % align);
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::MemoryPool()
throw()
{
  currentBlock_ = 0;
  currentSlot_ = 0;
  lastSlot_ = 0;
  freeSlots_ = 0;
  alloc_block_times = 0;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::MemoryPool(const MemoryPool& memoryPool)
throw()
{
  MemoryPool();
  alloc_block_times = 0;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
template<class U>
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::MemoryPool(const MemoryPool<U>& memoryPool)
throw()
{
  MemoryPool();
  alloc_block_times = 0;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::~MemoryPool()
throw()
{
  slot_pointer_ curr = currentBlock_;
  while (curr != 0) {
    slot_pointer_ prev = curr->next;
    operator delete(reinterpret_cast<void*>(curr));
    curr = prev;
  }
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::pointer
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::address(reference x)
const throw()
{
  return &x;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::const_pointer
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::address(const_reference x)
const throw()
{
  return &x;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
void
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::allocateBlock()
{
  // Allocate space for the new block and store a pointer to the previous one
  data_pointer_ newBlock = reinterpret_cast<data_pointer_>
                           (operator new(BlockSize));
  reinterpret_cast<slot_pointer_>(newBlock)->next = currentBlock_;
  currentBlock_ = reinterpret_cast<slot_pointer_>(newBlock);
  // Pad block body to staisfy the alignment requirements for elements
  data_pointer_ body = newBlock + sizeof(slot_pointer_);
  size_type bodyPadding = padPointer(body, sizeof(slot_type_));
  currentSlot_ = reinterpret_cast<slot_pointer_>(body + bodyPadding);
  lastSlot_ = reinterpret_cast<slot_pointer_>
              (newBlock + BlockSize - sizeof(slot_type_) + 1);
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::pointer
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::allocate(size_type, const_pointer)
{
  if (freeSlots_ != 0) {
    pointer result = reinterpret_cast<pointer>(freeSlots_);
    freeSlots_ = freeSlots_->next;
    return result;
  }
  else {
    if (currentSlot_ >= lastSlot_) {
      alloc_block_times ++;
      //printf("allocing block\n");
      printf("alloc_block_times: %d\n", alloc_block_times);
      allocateBlock();
    }
    return reinterpret_cast<pointer>(currentSlot_++);
  }
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline void
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::deallocate(pointer p, size_type)
{
  if (p != 0) {
    reinterpret_cast<slot_pointer_>(p)->next = freeSlots_;
    freeSlots_ = reinterpret_cast<slot_pointer_>(p);
  }
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::size_type
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::max_size()
const throw()
{
  size_type maxBlocks = -1 / BlockSize;
  return (BlockSize - sizeof(data_pointer_)) / sizeof(slot_type_) * maxBlocks;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline void
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::construct(pointer p, const_reference val)
{
  new (p) value_type (val);
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline void
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::construct(pointer p, constructor_type1 c1, constructor_type2 c2, constructor_type3 c3)
{
  new (p) value_type (c1, c2, c3);
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline void
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::destroy(pointer p)
{
  p->~value_type();
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::pointer
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::newElement(const_reference val)
{
  pointer result = allocate();
  construct(result, val);
  return result;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline typename MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::pointer
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::newElement(constructor_type1 c1, constructor_type2 c2, constructor_type3 c3)
{
  pointer result = allocate();
  construct(result, c1, c2, c3);
  return result;
}



template <typename T, size_t BlockSize, typename constructor_type1, typename constructor_type2, typename constructor_type3>
inline void
MemoryPool<T, BlockSize, constructor_type1, constructor_type2, constructor_type3>::deleteElement(pointer p)
{
  if (p != 0) {
    p->~value_type();
    deallocate(p);
  }
}



#endif // MEMORY_BLOCK_TCC
