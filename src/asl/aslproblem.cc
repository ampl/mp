/*
 A C++ interface to an AMPL problem.

 Copyright (C) 2012 AMPL Optimization Inc

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

#include "aslproblem.h"

#include <stdlib.h>
#include <cstring>

#ifndef _WIN32
# include <unistd.h>
#else
# include <io.h>
#endif

#include "mp/nl.h"
#include "mp/os.h"
#include "aslbuilder.h"
#include "expr-writer.h"

#ifdef _WIN32
# define close _close
#endif

#ifndef MP_HAVE_MKSTEMPS
extern "C" int mkstemps(char *pattern, int suffix_len);
#endif

using mp::asl::internal::ASLBuilder;

namespace {

// An .nl handler that builds an ASL problem using ASLBuilder.
class ASLHandler : public mp::ProblemBuilderToNLAdapter<ASLBuilder> {
 private:
  int flags_;

  typedef mp::ProblemBuilderToNLAdapter<ASLBuilder> Base;

 public:
  explicit ASLHandler(ASLBuilder &b, int obj_index = 0) : Base(b, obj_index) {}

  int flags() const { return flags_; }

  void OnHeader(const mp::NLHeader &h) {
    Base::OnHeader(h);
    flags_ = h.flags;
  }
};
}

namespace mp {

Solution::Solution()
: solve_code_(-1), num_vars_(0), num_cons_(0), values_(0), dual_values_(0) {}

Solution::~Solution() {
  std::free(values_);
  std::free(dual_values_);
}

void Solution::Swap(Solution &other) {
  std::swap(solve_code_, other.solve_code_);
  std::swap(num_vars_, other.num_vars_);
  std::swap(num_cons_, other.num_cons_);
  std::swap(values_, other.values_);
  std::swap(dual_values_, other.dual_values_);
}

void Solution::Read(fmt::StringRef stub, int num_vars, int num_cons) {
  // Allocate filename large enough to hold stub, ".sol" and terminating zero.
  std::vector<char> filename(stub.size() + 5);
  std::strcpy(&filename[0], stub.c_str());
  ASL asl = ASL();
  asl.i.n_var_ = num_vars;
  asl.i.n_con_ = num_cons;
  asl.i.ASLtype = 1;
  asl.i.filename_ = &filename[0];
  asl.i.stub_end_ = asl.i.filename_ + stub.size();
  Solution sol;
  sol.num_vars_ = num_vars;
  sol.num_cons_ = num_cons;
  char *message = read_sol_ASL(&asl, &sol.values_, &sol.dual_values_);
  if (!message)
    throw Error("Error reading solution file");
  std::free(message);
  Swap(sol);
  solve_code_ = asl.p.solve_code_;
}

void ASLProblem::Free() {
  if (var_capacity_) {
    delete [] asl_->i.LUv_;
    delete [] asl_->i.Uvx_;
    delete [] var_types_;
    asl_->i.LUv_ = asl_->i.Uvx_ = 0;
    var_capacity_ = 0;
  }
  ASL_fg *fg = asl_->i.ASLtype == ASL_read_fg ?
        reinterpret_cast<ASL_fg*>(asl_) : 0;
  if (obj_capacity_) {
    if (fg) {
      delete [] fg->I.obj_de_;
      fg->I.obj_de_ = 0;
    }
    delete [] asl_->i.objtype_;
    delete [] asl_->i.Ograd_;
    asl_->i.objtype_ = 0;
    asl_->i.Ograd_ = 0;
    obj_capacity_ = 0;
  }
  if (logical_con_capacity_) {
    if (fg) {
      delete [] fg->I.lcon_de_;
      fg->I.lcon_de_ = 0;
    }
    logical_con_capacity_ = 0;
  }
}

ASLProblem::ASLProblem()
: asl_(ASL_alloc(ASL_read_fg)), var_capacity_(0), obj_capacity_(0),
  logical_con_capacity_(0), var_types_(0) {
}

ASLProblem::ASLProblem(Proxy proxy)
: asl_(proxy.asl_), var_capacity_(0), obj_capacity_(0),
  logical_con_capacity_(0), var_types_(0) {
  proxy.asl_ = 0;
}

ASLProblem::~ASLProblem() {
  Free();
  ASL_free(reinterpret_cast<ASL**>(&asl_));
}

// A manager of temporary files.
class TempFiles {
 private:
  fmt::internal::MemoryBuffer<char, fmt::internal::INLINE_BUFFER_SIZE> name_;

  FMT_DISALLOW_COPY_AND_ASSIGN(TempFiles);

 public:
  TempFiles();

  ~TempFiles() {
    const char *stub = this->stub();
    std::remove(fmt::format("{}.nl", stub).c_str());
    std::remove(fmt::format("{}.sol", stub).c_str());
  }

  const char *stub() const { return &name_[0]; }
};

TempFiles::TempFiles() {
  std::string temp_dir = path::temp_directory_path().string();
  const char *s = &temp_dir[0];
  name_.append(s, s + temp_dir.size());
  const char TEMPLATE[] = "/XXXXXX.nl";
  name_.append(TEMPLATE, TEMPLATE + sizeof(TEMPLATE));  // include nul char
  const int SUFFIX_LEN = 3;  // length of the ".nl" suffix
  int fd = mkstemps(&name_[0], SUFFIX_LEN);
  if (fd == -1)
    throw fmt::SystemError(errno, "cannot create temporary file {}", &name_[0]);
  name_[name_.size() - SUFFIX_LEN] = '\0';  // remove .nl suffix
  close(fd);
}

ASLSuffixPtr SuffixView::Find(const char *name, unsigned flags) const {
  for (SufDesc *d = asl_->i.suffixes[kind_]; d; d = d->next) {
    if (!std::strcmp(name, d->sufname)) {
      return ASLSuffixPtr(asl_, (flags & suf::INPUT) != 0 &&
          (d->kind & suf::INPUT) == 0 ? 0 : d);
    }
  }
  return ASLSuffixPtr();
}

void ASLProblem::AddVar(double lb, double ub, var::Type type) {
  int &num_vars = asl_->i.n_var_;
  if (num_vars >= var_capacity_) {
    IncreaseCapacity(num_vars, var_capacity_);
    Grow(asl_->i.LUv_, num_vars, var_capacity_);
    Grow(asl_->i.Uvx_, num_vars, var_capacity_);
    if (var_types_)
      Grow(var_types_, num_vars, var_capacity_);
  }
  if (type != var::CONTINUOUS) {
    // Allocate var_types_ if this is the first integer variable added
    // after continuous.
    int num_integer_vars = ASLProblem::num_integer_vars();
    if (!var_types_ && num_vars != num_integer_vars) {
      var_types_ = new var::Type[var_capacity_];
      std::fill_n(fmt::internal::make_ptr(var_types_, var_capacity_),
        num_integer_vars, var::INTEGER);
      std::fill(var_types_ + num_integer_vars,
          var_types_ + num_vars, var::CONTINUOUS);
    }
    ++asl_->i.niv_;
  }
  asl_->i.LUv_[num_vars] = lb;
  asl_->i.Uvx_[num_vars] = ub;
  if (var_types_)
    var_types_[num_vars] = type;
  ++num_vars;
}

void ASLProblem::AddObj(obj::Type type, asl::NumericExpr expr) {
  int &num_objs = asl_->i.n_obj_;
  ASL_fg *fg = asl_->i.ASLtype == ASL_read_fg ?
        reinterpret_cast<ASL_fg*>(asl_) : 0;
  if (num_objs >= obj_capacity_) {
    IncreaseCapacity(num_objs, obj_capacity_);
    if (fg)
      Grow(fg->I.obj_de_, num_objs, obj_capacity_);
    Grow(asl_->i.objtype_, num_objs, obj_capacity_);
    Grow(asl_->i.Ograd_, num_objs, obj_capacity_);
  }
  cde e = cde();
  e.e = expr.impl_;
  if (fg)
    fg->I.obj_de_[num_objs] = e;
  asl_->i.objtype_[num_objs] = type;
  asl_->i.Ograd_[num_objs] = 0;
  ++num_objs;
  ++asl_->i.nlo_;
}

void ASLProblem::AddCon(asl::LogicalExpr expr) {
  if (asl_->i.ASLtype != ASL_read_fg)
    throw Error("problem doesn't support logical constraints");
  ASL_fg *fg = reinterpret_cast<ASL_fg*>(asl_);
  int &num_logical_cons = asl_->i.n_lcon_;
  if (num_logical_cons >= logical_con_capacity_) {
    IncreaseCapacity(num_logical_cons, logical_con_capacity_);
    Grow(fg->I.lcon_de_, num_logical_cons, logical_con_capacity_);
  }
  cde e = cde();
  e.e = expr.impl_;
  fg->I.lcon_de_[num_logical_cons] = e;
  ++num_logical_cons;
}

void ASLProblem::Read(fmt::StringRef stub, unsigned flags) {
  Free();

  // Add the .nl extension if necessary.
  const char EXT[] = ".nl";
  std::size_t ext_size = sizeof(EXT) - 1;
  fmt::MemoryWriter name;
  name << stub;
  if (name.size() < ext_size ||
      std::strcmp(name.c_str() + name.size() - ext_size, EXT) != 0) {
    name << EXT;
  }

  asl_->p.want_derivs_ = 0;
  asl_->i.want_xpi0_ = (flags & READ_INITIAL_VALUES) != 0;
  asl::internal::ASLBuilder builder(reinterpret_cast<ASL*>(asl_));
  builder.set_flags(ASL_allow_CLP | ASL_sep_U_arrays |
                    ASL_allow_missing_funcs |
                    asl::internal::ASL_STANDARD_OPCODES | flags);
  builder.set_stub(name.c_str());
  ASLHandler handler(builder, ASLHandler::NEED_ALL_OBJS);
  ReadNLFile(name.c_str(), handler);
  asl_->i.flags = handler.flags();
}

void ASLProblem::WriteNL(
    fmt::StringRef stub, ProblemChanges *pc, unsigned flags) {
  int nfunc = asl_->i.nfunc_;
  if ((flags & IGNORE_FUNCTIONS) != 0)
    asl_->i.nfunc_ = 0;
  int result = fg_write_ASL(reinterpret_cast<ASL*>(asl_),
      stub.c_str(), pc ? pc->vco() : 0, ASL_write_ASCII);
  asl_->i.nfunc_ = nfunc;
  if (result)
    throw Error("Error writing .nl file");
}

void ASLProblem::Solve(fmt::StringRef solver_name,
    Solution &sol, ProblemChanges *pc, unsigned flags) {
  TempFiles temp;
  WriteNL(temp.stub(), pc, flags);
  // Run the solver and read the solution file.
  fmt::MemoryWriter command;
  command.write("{} {} -AMPL", solver_name.c_str(), temp.stub());
  int exit_code = std::system(command.c_str());
  if (exit_code != 0) {
    throw Error(
        "Error running solver {}, exit code = {}", solver_name, exit_code);
  }
  sol.Read(temp.stub(), num_vars() + (pc ? pc->num_vars() : 0),
      num_cons() + (pc ? pc->num_cons() : 0));
}

NewVCO *ProblemChanges::vco() {
  static double dummy;
  vco_.nnv = static_cast<int>(var_lb_.size());
  vco_.nnc = static_cast<int>(cons_.size());
  vco_.nno = static_cast<int>(objs_.size());
  if (!var_lb_.empty()) {
    vco_.LUnv = &var_lb_[0];
    vco_.Unv = &var_ub_[0];
  } else {
    vco_.LUnv = &dummy;
  }
  if (!cons_.empty()) {
    vco_.LUnc = &con_lb_[0];
    vco_.Unc = &con_ub_[0];
    vco_.newc = &cons_[0];
  } else {
    vco_.LUnc = &dummy;
  }
  if (!objs_.empty()) {
    vco_.newo = &objs_[0];
    vco_.ot = &obj_types_[0];
  }
  return &vco_;
}

void ProblemChanges::AddObj(
    obj::Type type, unsigned size, const double *coefs, const int *vars) {
  std::size_t start = obj_terms_.size();
  obj_terms_.resize(start + size);
  ograd dummy;
  ograd *prev = &dummy;
  for (unsigned i = 0; i < size; ++i) {
    ograd &term = obj_terms_[start + i];
    term.coef = coefs[i];
    term.varno = vars[i];
    prev->next = &term;
    prev = &term;
  }
  objs_.push_back(&obj_terms_[start]);
  obj_types_.push_back(type);
}

void ProblemChanges::AddCon(const double *coefs, double lb, double ub) {
  con_lb_.push_back(lb);
  con_ub_.push_back(ub);
  std::size_t start = con_terms_.size();
  int num_vars = static_cast<int>(problem_->num_vars() + var_lb_.size());
  con_terms_.resize(start + num_vars);
  ograd dummy;
  ograd *prev = &dummy;
  for (int i = 0; i < num_vars; ++i) {
    ograd &term = con_terms_[start + i];
    term.coef = coefs[i];
    term.varno = i;
    prev->next = &term;
    prev = &term;
  }
  cons_.push_back(&con_terms_[start]);
}

ProblemChanges::ProblemChanges(const ProblemChanges &other) {
  *this = other;
}

ProblemChanges &ProblemChanges::operator=(const ProblemChanges &rhs) {
  if (this == &rhs) {
    return *this;
  }
  problem_ = rhs.problem_;
  var_lb_ = rhs.var_lb_;
  var_ub_ = rhs.var_ub_;
  con_lb_ = rhs.con_lb_;
  con_ub_ = rhs.con_ub_;
  con_terms_ = rhs.con_terms_;
  obj_terms_ = rhs.obj_terms_;
  obj_types_ = rhs.obj_types_;
  cons_.resize(rhs.cons_.size());
  objs_.resize(rhs.objs_.size());
  vco_ = rhs.vco_;

  int next = 0;
  for (size_t i = 0; i < rhs.con_terms_.size(); ++i) {
    if (next < rhs.num_cons() && rhs.cons_[next] == &(rhs.con_terms_[i])) {
      cons_[next] = &(con_terms_[i]);
      ++next;
    }
    if (con_terms_[i].next) {
      con_terms_[i].next = &(con_terms_[i+1]);
    }
  }

  next = 0;
  for (size_t i = 0; i < rhs.obj_terms_.size(); ++i) {
    if (next < rhs.num_objs() && rhs.objs_[next] == &(rhs.obj_terms_[i])) {
      objs_[next] = &(obj_terms_[i]);
      ++next;
    }
    if (obj_terms_[i].next) {
      obj_terms_[i].next = &(obj_terms_[i+1]);
    }
  }
  return *this;
}

fmt::Writer &operator<<(fmt::Writer &w, const ASLProblem &p) {
  Write(w, p);
}
}
