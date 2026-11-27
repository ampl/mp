var x{1..3} integer >= -2 <= 1;

maximize o: x[1] + 4*x[2] + 5*x[3];

s.t. c1: numberof -2 in (x[1], x[2], x[3]) == 0;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: numberof -2 in (x[1], x[2], x[3]) == 0;
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
