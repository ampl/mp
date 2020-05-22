#ifndef INTERFACEMIPTEST_H
#define INTERFACEMIPTEST_H

#include <vector>

#include "mp/convert/MIP/mp2mip.h"

namespace interface_test {

enum Sense {
  minimize_ = mp::obj::MIN,
  maximize_ = mp::obj::MAX
};
enum VarType {
  I_ = mp::var::INTEGER,
  F_ = mp::var::CONTINUOUS
};
constexpr double infty_ = 1e100;

struct MIPInstance {
  struct SparseVec {       // TODO Use LinearExpr for this
    std::vector<double> c_;
    std::vector<int>    v_;
    SparseVec(const std::vector<double>& c, const std::vector<int>& v) :
      c_(c), v_(v) { }
    SparseVec(int nnz, const double* c, const int* v) :
      c_(c, c+nnz), v_(v, v+nnz) { }
    bool operator==(const SparseVec& sv) const { return c_==sv.c_ && v_==sv.v_; }
    size_t size() const { return c_.size(); }
  };
  struct Objective {
    Sense sense_;
    SparseVec le_;    // lin expression. TODO non-lin as well
  };
  using ObjectivesContainer = std::vector<Objective>;
  ObjectivesContainer objs_;
  /// Variables
  std::vector<double> varLBs_, varUBs_;
  std::vector<VarType> varTypes_;

  struct Constraint {
    SparseVec le_;
    double lb_;
    double ub_;
  };
  using ConstraintsContainer = std::vector<Constraint>;
  ConstraintsContainer cons_;

  bool ApproxEqual(const MIPInstance& mip) const {
    return
        objs_.size() == mip.objs_.size() &&
        varLBs_.size() == mip.varLBs_.size() &&
        cons_.size() == mip.cons_.size();
  }
};

template <class Interface>
void feedInstance( Interface& interface, const MIPInstance& mip ) {
  interface.InputVariables(mip.varLBs_.size(),
                           mip.varLBs_.data(), mip.varUBs_.data(),
                           (const mp::var::Type*)mip.varTypes_.data());
  for (MIPInstance::ObjectivesContainer::const_iterator it = mip.objs_.begin();
       it!=mip.objs_.end(); ++it) {
    interface.InputObjective((mp::obj::Type)it->sense_,
                             it->le_.size(), it->le_.c_.data(), it->le_.v_.data());
  }
  for (MIPInstance::ConstraintsContainer::const_iterator it = mip.cons_.begin();
       it!=mip.cons_.end(); ++it) {
    interface.InputAlgebraicCon(it->le_.size(), it->le_.c_.data(), it->le_.v_.data(),
                             it->lb_, it->ub_);
  }
}

/// A toy backend using struct MIPInstance
class MIPInstanceBackend :
    public mp::BasicBackend<MIPInstanceBackend>
{
  MIPInstance instance_;
public:
  MIPInstanceBackend() : mp::BasicBackend<MIPInstanceBackend>("MIPInstanceBackend") { }
  MIPInstance& GetInstance() { return instance_; }

  /// These things the concrete interface currently has to define
  void AddVariables(int n, double* lbs, double* ubs, mp::var::Type* types) {
    instance_.varLBs_.insert(instance_.varLBs_.end(), lbs, lbs+n);
    instance_.varUBs_.insert(instance_.varUBs_.end(), ubs, ubs+n);
    instance_.varTypes_.insert(instance_.varTypes_.end(),
                               (VarType*)types, (VarType*)types+n);
  }
  void AddLinearObjective(mp::obj::Type sense, int nnz,
                          const double* c, const int* v) {
    instance_.objs_.push_back({(Sense)sense, { nnz, c, v }});
  }

  void AddLinearConstraint(int nnz, const double* c, const int* v,
                           double lb, double ub) {
    instance_.cons_.push_back({ { nnz, c, v }, lb, ub });
  }

};

/// Testing the default MIP interface layer
class MIPInterfaceTester :
    public mp::MPToMIPConverter<MIPInterfaceTester, MIPInstanceBackend>
{
public:
  /// This is testing API
  bool OutputModelSeemsEqualTo(const MIPInstance& mip) {
    return GetBackend().GetInstance().ApproxEqual( mip );
  }
};

} // namespace interface_test


#endif // INTERFACEMIPTEST_H
