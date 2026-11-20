var x{1..3} integer >= 0 <= 2;

minimize o: 2*x[1] - 5*x[2] - x[3];

s.t. c1: not ((x[2] != 0) and (x[1] == 2));
s.t. c2: abs(x[1] - x[1]) > 2;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
