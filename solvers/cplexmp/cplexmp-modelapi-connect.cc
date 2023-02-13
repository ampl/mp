#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

#include "cplexmpmodelapi.h"


namespace mp {

/// Defining the function in ...modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateCplexModelMgr(CplexCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			CplexModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
