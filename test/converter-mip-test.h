#ifndef CONVERTERMIPTEST_H
#define CONVERTERMIPTEST_H

#include "mp/arrayref.h"

#include "mp/flat/model_flattener.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/backend_model_api_base.h"

namespace mip_converter_test {

enum Sense {
  minimize_ = mp::obj::MIN,
  maximize_ = mp::obj::MAX
};
constexpr auto I_ = mp::var::INTEGER;
constexpr auto F_ = mp::var::CONTINUOUS;
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
    SparseVec le_;        // lin expression
    bool operator==(const Objective& obj) const {
      bool eq1 = sense_==obj.sense_;
      bool eq2 = le_==obj.le_;
      if (!eq1)
        printf("Diff sense\n");
      if (!eq2) {
        printf("Diff coefs\n");
        for (size_t i=0; i<le_.size(); ++i) {
          printf(" %d:%f", le_.v_[i], le_.c_[i]);
        }
        printf("\n");
        for (size_t i=0; i<obj.le_.size(); ++i) {
          printf(" %d:%f", obj.le_.v_[i], obj.le_.c_[i]);
        }
        printf("\n");
      }
      return eq1 && eq2;
    }
  };
  using ObjectivesContainer = std::vector<Objective>;
  ObjectivesContainer objs_;
  /// Variables
  mp::VarArrayDef vars_;

  struct Constraint {
    SparseVec le_;
    double lb_;
    double ub_;
  };
  using ConstraintsContainer = std::vector<Constraint>;
  ConstraintsContainer cons_;

  bool ObjsEqual(const MIPInstance& mip) const {
    return
        objs_ == mip.objs_;
  }
  bool VarBoundsEqual(const MIPInstance& mip) const {
    return
        vars_.size() == mip.vars_.size() &&
        std::equal(vars_.plb(), vars_.plb()+vars_.size(),
                   mip.vars_.plb()) &&
        std::equal(vars_.pub(), vars_.pub()+vars_.size(),
                   mip.vars_.pub());
  }
  bool NConstrEqual(const MIPInstance& mip) const {
    return
        cons_.size() == mip.cons_.size();
  }
};

template <class Interface>
void feedInstance( Interface& interface, const MIPInstance& mip ) {
  interface.InputVariables(mip.vars_.size(),
                           mip.vars_.plb(), mip.vars_.pub(),
                           mip.vars_.ptype());
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
    public mp::BasicBackendFlatModelAPI,
    public mp::EnvKeeper
{
  MIPInstance instance_;
public:
  MIPInstanceBackend(mp::Env& e) : mp::EnvKeeper(e) { }

  static constexpr const char* GetName() { return "tester"; }

  MIPInstance& GetInstance() { return instance_; }

  /// These things the concrete interface currently has to define
  void AddVariables(const mp::VarArrayDef& v) {
    instance_.vars_ = v;
  }
  void SetLinearObjective(int , const mp::LinearObjective& lo) {
    mip_converter_test::MIPInstance::SparseVec lin_part {lo.coefs(), lo.vars()};
    instance_.objs_.push_back({(Sense)lo.obj_sense(),
                              std::move(lin_part)});
  }

  /// Allow all constraint types to be compiled
  USE_BASE_CONSTRAINT_HANDLERS(mp::BasicBackendFlatModelAPI)

  ACCEPT_CONSTRAINT(mp::LinConEQ, mp::Recommended, mp::CG_Default)

  /// Specialize for LinearConstraint
  void AddConstraint(const mp::LinConEQ& lc) {
    instance_.cons_.push_back({ { (int)lc.size(), lc.pcoefs(), lc.pvars() },
                                lc.lb(), lc.ub() });
  }
};

/// Testing the default MIP interface layer
class MIPConverterTester :
    public mp::ModelFltImpl<mp::ModelFlattener, mp::Problem,
      mp::FlatCvtImpl<mp::MIPFlatConverter, MIPInstanceBackend> >
//    public mp::MPToMIPConverter<MIPConverterTester, MIPInstanceBackend>
{
  using Base = mp::ModelFltImpl<mp::ModelFlattener, mp::Problem,
    mp::FlatCvtImpl<mp::MIPFlatConverter, MIPInstanceBackend> >;
public:
  /// Construct
  MIPConverterTester(mp::Env& e) : Base(e) { }
  /// This is testing API
  bool ObjsEqual(const MIPInstance& mip) {
    return GetBackend().GetInstance().ObjsEqual( mip );
  }
  bool VarBoundsEqual(const MIPInstance& mip) {
    return GetBackend().GetInstance().VarBoundsEqual( mip );
  }
  bool NConstrEqual(const MIPInstance& mip) {
    return GetBackend().GetInstance().NConstrEqual( mip );
  }
protected:
  MIPInstanceBackend& GetBackend() {
    return GetFlatCvt().GetBasicBackend();
  }
};

} // namespace mip_converter_test


#endif // CONVERTERMIPTEST_H
