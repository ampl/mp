#ifndef FLAT_MODEL_API_CONNECT_H
#define FLAT_MODEL_API_CONNECT_H

#include <memory>

#include "mp/env.h"
#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/converter.h"
#include "mp/flat/problem_flattener.h"

namespace mp {

/// A template to create a ModelManager with a FlatConverter
/// for a specific solver.
/// Currently is uses mp::Problem as intermediate storage
/// of the NL model.
///
/// @tparam Backend: the Backend class
/// @tparam ModelAPI: the ModelAPI class to be created and
///   used by the ModelManager
/// @tparam FlatConverter: the FlatConverter to be used
///
/// @param gc: Backend object
/// @param e: the MP environment
/// @param pPre: return ValuePresolver for Backends that use it
template <class ModelAPI,
          template <typename, typename, typename> class FlatConverter,
          class Backend>
std::unique_ptr<BasicModelManager>
CreateModelMgrWithFlatConverter(Backend& gc, Env& e,
                     pre::BasicValuePresolver*& pPre) {
  using SolverFlatCvt = FlatCvtImpl<FlatConverter, ModelAPI>;
  using SolverProblemFlattener = mp::ProblemFltImpl<
    mp::ProblemFlattener, mp::Problem, SolverFlatCvt>;
  auto pcvt = new SolverProblemFlattener(e);
  auto res = CreateModelManagerWithStdBuilder(
        std::unique_ptr< BasicConverter<mp::Problem> >{ pcvt } );
  pcvt->GetFlatCvt().GetModelAPI().set_other(&gc);
  gc.set_other(
        &pcvt->GetFlatCvt().GetModelAPI());
  pPre = &pcvt->GetFlatCvt().GetValuePresolver();
  return res;
}

} // namespace mp

#endif // FLAT_MODEL_API_CONNECT_H
