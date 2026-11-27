var x{1..3} integer >= 0 <= 3;

maximize o: x[1] - 2*x[2] - 4*x[3];

s.t. c1: numberof 4 in (x[1], x[2], x[3]) <= 1;
s.t. c2: exists {i in 1..3} x[i] == 3;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: ((x[2] > -1) ==> (x[1] > 3)) and ((x[1] > 3) ==> (x[1] >= -1));
