######################################################
## names_const_01.mod: shared constant variables are
## named when the model passes names (auxfiles rc).
## The indicator's conclusion x >= 2 is false by the
## bounds of x, and the if-then-else has the implicit
## else 0: MP represents both with constant variables.
######################################################

var x binary;
var u binary;
var z >= 0, <= 10;

maximize O: x + u + z;

s.t. Ind: u >= 0.5 ==> x >= 2;
s.t. Ite: z <= if x >= 1 then 5;
