#include "highsmpmodelapi.h"

#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in.modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateHighsModelMgr(HighsCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			HighsModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
