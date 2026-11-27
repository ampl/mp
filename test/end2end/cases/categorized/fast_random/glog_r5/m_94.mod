var x{1..3} integer >= 0 <= 4;

maximize o: 5*x[1] + 5*x[2] - x[3];

s.t. c1: numberof 4 in (x[1], x[2], x[3]) == 0;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: numberof 2 in (x[1], x[2], x[3]) > 2;
