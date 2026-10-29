#ifndef CBCMPMODELAPI_H
#define CBCMPMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "cbcmpcommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"


namespace mp {

class CbcmpModelAPI :
    public CbcmpCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

  // Struct to store the model being generated
  struct ModelInstance {

    // Store the row representation of the model and offers conversion to CSC 
    class SparseMatrixBuilder {
      int num_cols;
      std::vector<std::vector<std::pair<int, double>>> rows; 

    public:
      void addrow(const int* varindices, const double* coefficients, int size) {
        std::vector<std::pair<int, double>> row;
        for (int i = 0; i < size; ++i) {
          row.emplace_back(varindices[i], coefficients[i]);
        }
        rows.push_back(row);
      }

      void toCSC(std::vector<double>& values, std::vector<int>& row_indices, std::vector<CoinBigIndex>& col_pointers) {
        std::vector<std::vector<std::pair<int, double>>> cols(num_cols);

        // Transpose rows to columns
        for (int row = 0; row < rows.size(); ++row) {
          for (const auto& entry : rows[row]) {
            int col = entry.first;
            double value = entry.second;
            cols[col].emplace_back(row, value);
          }
        }

        // Populate CSC arrays
        col_pointers.push_back(0);
        for (int col = 0; col < num_cols; ++col) {
          for (const auto& entry : cols[col]) {
            row_indices.push_back(entry.first);
            values.push_back(entry.second);
          }
          col_pointers.push_back(values.size());
        }
      }
      void setNumCols(int nc) { num_cols = nc; }
    };

    struct SOS {
      std::vector<int> sos_rowstarts;
      std::vector<int> sos_colindices;
      std::vector<double> sos_colweights;
      std::vector<int> sos_types;
      void Clear() {
        sos_rowstarts = std::vector<int>();
        sos_colindices = std::vector<int>();
        sos_colweights = std::vector<double>();
        sos_types = std::vector<int>();
      }
    };

    // All model info
    SparseMatrixBuilder matrix;
    int objsense;
    std::vector<std::string> varnames, connames;
    std::vector<double> lb, ub, lhs, rhs;
    std::vector<double> obj;
    std::vector<int> varinteger;
    SOS soslist[2];

    bool visited = false;

    void addSOS(std::size_t size, const int* indices,
      const double* weights, int type) {
      soslist[type-1].sos_rowstarts.push_back(soslist[type - 1].sos_colindices.size());
      soslist[type - 1].sos_colindices.insert(soslist[type - 1].sos_colindices.end(), indices, indices + size);
      soslist[type - 1].sos_colweights.insert(soslist[type - 1].sos_colweights.end(), weights, weights + size);
      soslist[type - 1].sos_types.push_back(type);
    }

    void SetNumVars(int n) {
      varnames.resize(n);
      lb.resize(n);
      ub.resize(n);
      obj.resize(n);
      matrix.setNumCols(n);
    }

    void AddRow(const char* name, std::size_t size, const int* vars,
      const double* coefs, double lhss, double rhss) {
      connames.push_back(name);
      matrix.addrow(vars, coefs, size);
      lhs.push_back(lhss);
      rhs.push_back(rhss);
    }
    void Clear() {
      varnames = std::vector<std::string>();
      connames = std::vector<std::string>();
      lb = std::vector<double>();
      ub = std::vector<double>();
      lhs = std::vector<double>();
      rhs = std::vector<double>();
      obj = std::vector<double>();
      varinteger = std::vector<int>();
      for (int i = 0; i < 2; i++)
        soslist[i].Clear();

    }
  };

  ModelInstance m;
public:
  /// Constructor
  CbcmpModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CbcmpModelAPI"; }

  /// Called before problem input.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After problem input.
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }


  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// TODO For each suppoted constraint type, add the ACCEPT_CONSTRAINT macro
  /// and the relative AddConstraint function.
  /// Below some typical constraint handlers of a MIP solver.
  /// Further constraint types which could be handled natively by some solvers:
  /// - IndicatorConstraint(Lin/Quad)(LE/EQ/GE)
  /// - Multidirectional indicators Cond(Lin/Quad)Con(LT/LE/EQ/GE/GT), where
  ///   the implication direction (</==/>) depends in the context
  /// - Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constr_std.h and other drivers.


  /// The linear range constraint, if fully supported with basis info etc.
  //ACCEPT_CONSTRAINT(LinConRange, NotAccepted, CG_Linear)

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
    void AddConstraint(const LinConRange& lc);
  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// Set ``option pl_linearize 0;`` in AMPL if the solver
  /// supports PL natively.
  /// Not recommended: CBC 2.10.12 returns wrong solutions with some
  /// SOS models; MP linearizes them by default. acc:sos1/sos2=2
  /// restore passing them to CBC.
  ACCEPT_CONSTRAINT(SOS1Constraint, AcceptedButNotRecommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, AcceptedButNotRecommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

};

} // namespace mp

#endif // CBCMPMODELAPI_H
