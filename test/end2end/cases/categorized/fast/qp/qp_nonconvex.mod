##

var x >= 0;
var z >= 0;
var y >= 1 integer;

maximize DivRedefined: x*z;

s.t. RedefDiv: z*y <= 1;

s.t. Lin15: x+y+z <= 19;
