var x{1..3} integer >= 0 <= 3;

minimize o: -4*x[1] - 5*x[2] + 0*x[3];

s.t. c1: exists {i in 1..3} x[i] == 3;
s.t. c2: (x[2] < 0) ==> (x[2] > 3);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
