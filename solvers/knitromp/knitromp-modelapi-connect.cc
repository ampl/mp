#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

#include "knitrompmodelapi.h"


namespace mp {

/// Defining the function in ...-modelapi-connect.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateKnitrompModelMgr(KnitrompCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			KnitrompModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
