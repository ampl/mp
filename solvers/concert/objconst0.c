
#define SKIP_NL2_DEFINES
#include "nlp.h"
#undef f_OPNUM
#include "r_opn0.hd"

#ifdef __cplusplus
extern "C"
#endif
 real
objconst0(ASL_fg *a, int n)
{
	expr_n *e = (expr_n*)(a->I.obj_de_ + n)->e;
	if (e->op == (efunc_n*)(unsigned long)f_OPNUM)
		return e->v;
	return 0;
	}
