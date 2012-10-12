#include "gecode.h"

#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>

#include <memory>

using namespace Gecode;

class GecodeModel : public Space {
 private:
  IntVarArray l;
 public:
  GecodeModel() : l(*this, 8, 0, 9) {
    // TODO: convert from AMPL to Gecode format
    IntVar s(l[0]), e(l[1]), n(l[2]), d(l[3]),
           m(l[4]), o(l[5]), r(l[6]), y(l[7]);
    // no leading zeros
    rel(*this, s, IRT_NQ, 0);
    rel(*this, m, IRT_NQ, 0);
    // all letters distinct
    distinct(*this, l);
    // linear equation
    IntArgs c(4+4+5); IntVarArgs x(4+4+5);
    c[0]=1000; c[1]=100; c[2]=10; c[3]=1;
    x[0]=s;    x[1]=e;   x[2]=n;  x[3]=d;
    c[4]=1000; c[5]=100; c[6]=10; c[7]=1;
    x[4]=m;    x[5]=o;   x[6]=r;  x[7]=e;
    c[8]=-10000; c[9]=-1000; c[10]=-100; c[11]=-10; c[12]=-1;
    x[8]=m;      x[9]=o;     x[10]=n;    x[11]=e;   x[12]=y;
    linear(*this, c, x, IRT_EQ, 0);
    // post branching
    branch(*this, l, INT_VAR_SIZE_MIN, INT_VAL_MIN);
  }
  // search support
  GecodeModel(bool share, GecodeModel &s) : Space(share, s) {
    l.update(*this, share, s.l);
  }
  virtual Space *copy(bool share) {
    return new GecodeModel(share,*this);
  }
  // print solution
  void print(void) const {
    std::cout << l << std::endl;
  }
};

Driver::~Driver() {}

void Driver::run(char **argv) {
  // create model and search engine
  std::auto_ptr<GecodeModel> model(new GecodeModel);
  DFS<GecodeModel> e(model.get());
  model.reset();
  // search and print all solutions
  for (;;) {
    std::auto_ptr<GecodeModel> s(e.next());
    if (!s.get()) break;
    s->print();
  }
}
