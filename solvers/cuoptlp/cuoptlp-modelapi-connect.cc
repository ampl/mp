#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

#include "cuoptlpmodelapi.h"


namespace mp {

/// Defining the function in ...-modelapi-connect.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateCuoptlpModelMgr(CuoptlpCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			CuoptlpModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
