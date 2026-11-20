var x{1..3} integer >= 0 <= 4;

minimize o: 5*x[1] + 6*x[2] + 0*x[3];

s.t. c1: (x[2] != 4) <==> (x[2] + x[1] != 1);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: abs(x[2] - x[1]) != 0;
