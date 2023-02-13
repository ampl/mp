#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

#include "cbcmpmodelapi.h"


namespace mp {

std::unique_ptr<BasicModelManager>
CreateCbcmpModelMgr(CbcmpCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			CbcmpModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
