#ifndef MODEL_MGR_WITH_STD_PB_H
#define MODEL_MGR_WITH_STD_PB_H

/**
 * Generate ModelManagerWithPB<mp::Problem>, header file
 */
#include <memory>

#include "mp/model-mgr-base.h"
#include "mp/problem.h"
#include "mp/converter-base.h"

namespace mp {

/// Declare an ModelManager factory
std::unique_ptr<BasicModelManager>
CreateModelManagerWithStdBuilder(
    std::unique_ptr<BasicConverter<mp::Problem> > pcvt);

}

#endif // MODEL_MGR_WITH_STD_PB_H
