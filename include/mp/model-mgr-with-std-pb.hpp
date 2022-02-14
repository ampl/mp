#ifndef MODEL_MGR_WITH_STD_PB_HPP
#define MODEL_MGR_WITH_STD_PB_HPP

/**
 * Generate ModelManagerWithPB<mp::Problem>
 *
 * Include this in a separate .cc to improve compilation speed
 */
#include "mp/model-mgr-with-std-pb.h"
#include "mp/model-mgr-with-pb.h"

namespace mp {

using ModelMgrWithStdPB =
  ModelManagerWithProblemBuilder< BasicConverter< mp::Problem > >;

std::unique_ptr<BasicModelManager>
CreateModelManagerWithStdBuilder(
    std::unique_ptr<BasicConverter<mp::Problem> > pcvt) {
  return std::unique_ptr<BasicModelManager>
    { new ModelMgrWithStdPB(std::move(pcvt)) };
}

} // namespace mp

#endif // MODEL_MGR_WITH_STD_PB_HPP
