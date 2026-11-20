var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] - 2*x[2] + 2*x[3];

s.t. c1: numberof 0 in (x[1], x[2], x[3]) != 2;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: x[1] == 1;
s.t. c4: x[3] == 4;
