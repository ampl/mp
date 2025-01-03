#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

#include "cuoptmodelapi.h"


namespace mp {

/// Defining the function in ...-modelapi-connect.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateCuoptModelMgr(CuoptCommon& cc, Env& e,
										 pre::BasicValuePresolver*& pPre) {
	return CreateModelMgrWithFlatConverter<
			CuoptModelAPI, MIPFlatConverter >(cc, e, pPre);
}

}  // namespace mp
