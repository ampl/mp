/*
 Mock allocator

 Copyright (c) 2014, Victor Zverovich
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FMT_MOCK_ALLOCATOR_H_
#define FMT_MOCK_ALLOCATOR_H_

#include "gmock/gmock.h"

class MockAllocator {
 public:
  MockAllocator() {}
  MockAllocator(const MockAllocator &) {}

  MOCK_METHOD1(allocate, void* (std::size_t n));
  MOCK_METHOD2(deallocate, void (void* p, std::size_t n));
};

// A reference to a mock allocator which can be passed by value but point
// to the same underlying allocator used for testing.
template <typename Alloc = MockAllocator, typename T = char>
class AllocatorRef {
 private:
  Alloc *alloc_;

 public:
  typedef T value_type;

  template <typename U>
  struct rebind {
    typedef AllocatorRef<Alloc, U> other;
  };

  explicit AllocatorRef(Alloc *alloc = 0) : alloc_(alloc) {}

  AllocatorRef(const AllocatorRef &other) : alloc_(other.get()) {}

  template <typename U>
  AllocatorRef(const AllocatorRef<Alloc, U> &other) : alloc_(other.get()) {}

  AllocatorRef& operator=(const AllocatorRef &other) {
    alloc_ = other.alloc_;
    return *this;
  }

#if FMT_USE_RVALUE_REFERENCES
 private:
  void move(AllocatorRef &other) {
    alloc_ = other.alloc_;
    other.alloc_ = 0;
  }

 public:
  AllocatorRef(AllocatorRef &&other) {
    move(other);
  }

  AllocatorRef& operator=(AllocatorRef &&other) {
    assert(this != &other);
    move(other);
    return *this;
  }
#endif

  Alloc *get() const { return alloc_; }

  value_type* allocate(std::size_t n) {
    return reinterpret_cast<value_type*>(
          alloc_->allocate(n * sizeof(value_type)));
  }
  void deallocate(value_type* p, std::size_t n) { alloc_->deallocate(p, n); }
};

#endif  // FMT_MOCK_ALLOCATOR_H_
