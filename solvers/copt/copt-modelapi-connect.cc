#include "coptmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in ...-modelapi-connect.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateCoptModelMgr(CoptCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			CoptModelAPI, MIPFlatConverter >(cc, e, pPre);
}

} // namespace mp
