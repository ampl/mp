/*
 NL benchmarks

 Copyright (C) 2016 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "mp/nl-reader.h"
#include "benchmark/benchmark.h"

namespace obj = mp::obj;

struct Obj {
  obj::Type type;
  int expr;
  Obj(obj::Type t = obj::MIN, int e = 0) : type(t), expr(e) {}
};

// An NL handler that adds all objectives at once.
class SimultaneousObjHandler : public mp::NullNLHandler<int> {
 private:
  Obj *objs_;

 public:
  SimultaneousObjHandler() : objs_(0) {}
  ~SimultaneousObjHandler() { delete [] objs_; }

  void OnHeader(const mp::NLHeader &h) {
    objs_ = new Obj[h.num_objs]();
  }

  void OnObj(int index, obj::Type type, NumericExpr expr) {
    objs_[index] = Obj(type, expr);
  }
};

// An NL handler that adds objectives incrementally.
class IncrementalObjHandler : public mp::NullNLHandler<int> {
 private:
  std::vector<Obj> objs_;
  int obj_count_;

 public:
  IncrementalObjHandler() : obj_count_(0) {}

  void OnHeader(const mp::NLHeader &h) {
    objs_.reserve(h.num_objs);
  }

  void OnObj(int index, obj::Type type, NumericExpr expr) {
    if (index != obj_count_)
      throw std::runtime_error("invalid index");
    ++obj_count_;
    objs_.push_back(Obj(type, expr));
  }
};

// An NL handler that adds objectives incrementally using an array instead of
// a vector.
class IncrementalArrayObjHandler : public mp::NullNLHandler<int> {
 private:
  Obj *objs_;
  int obj_count_;

 public:
  IncrementalArrayObjHandler() : objs_(0), obj_count_(0) {}
  ~IncrementalArrayObjHandler() { delete [] objs_; }

  void OnHeader(const mp::NLHeader &h) {
    objs_ = new Obj[h.num_objs];
  }

  void OnObj(int index, obj::Type type, NumericExpr expr) {
    if (index != obj_count_)
      throw std::runtime_error("invalid index");
    objs_[obj_count_++] = Obj(type, expr);
  }
};

enum { NUM_OBJS = 10 };

static void AddObjsAtOnce(benchmark::State &state) {
  while (state.KeepRunning()) {
    SimultaneousObjHandler handler;
    mp::NLHeader header = mp::NLHeader();
    header.num_objs = NUM_OBJS;
    handler.OnHeader(header);
    for (int i = 0; i < NUM_OBJS; ++i)
      handler.OnObj(i, obj::MAX, 42);
  }
}

BENCHMARK(AddObjsAtOnce);

static void AddObjsIncVector(benchmark::State &state) {
  while (state.KeepRunning()) {
    IncrementalObjHandler handler;
    mp::NLHeader header = mp::NLHeader();
    header.num_objs = NUM_OBJS;
    handler.OnHeader(header);
    for (int i = 0; i < NUM_OBJS; ++i)
      handler.OnObj(i, obj::MAX, 42);
  }
}

BENCHMARK(AddObjsIncVector);

static void AddObjsIncArray(benchmark::State &state) {
  while (state.KeepRunning()) {
    IncrementalArrayObjHandler handler;
    mp::NLHeader header = mp::NLHeader();
    header.num_objs = NUM_OBJS;
    handler.OnHeader(header);
    for (int i = 0; i < NUM_OBJS; ++i)
      handler.OnObj(i, obj::MAX, 42);
  }
}

BENCHMARK(AddObjsIncArray);

BENCHMARK_MAIN()
