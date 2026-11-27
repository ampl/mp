# Exact integer equality must keep its RHS near a domain bound.
param k := 562949953421312;

var x integer >= k - 2, <= k + 1;
var b binary;

maximize Obj: b;

subject to Equiv: b <==> (x = k);
