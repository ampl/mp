 /**
  * Test expression map as well as if-then.
  * Test solution check and cap:eps.
  * Used in a Colab notebook on solution check.
  */

var x >=0, <= 1;
var y >=0, <= 1;
var b: binary;

s.t. ConImpl:
    b ==> 2*x + 3*y <= 5;

minimize TotalIf:
    if 2*x+3*y>5 then 2*x+3*y-3*b-25 else 2*x+3*y-3*b;