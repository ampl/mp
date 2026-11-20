# A small smooth model, solved by Knitro's Hessian options without exact
# second derivatives (hessopt 2, 3, 4, 6) with storequadcoefs=0: the
# driver used to stop with suf_get("xhonorbnd") fails.
var x {1..2} := -1.2;

minimize f: 100*(x[2] - x[1]^2)^2 + (1 - x[1])^2;

s.t. disk: x[1]^2 + x[2]^2 <= 2;
