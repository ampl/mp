#####################################################
# Test that the 2nd objective in the MO Emulator
# cleans up QP part
#####################################################

var x;
var y;

s.t. Con1: x+y <= 2;
s.t. Con1a: 3*x+y <= 5;
s.t. Con2: x-y <= 2;
s.t. Con2a: 3*x-y <= 5;

suffix objpriority IN;
suffix objabstol IN;

maximize Obj1QP: 2*x + y - 10*(x-2)^2
  suffix objpriority 5,
  suffix objabstol 1;
maximize Obj2Lin: 2*x - y suffix objpriority 2;
