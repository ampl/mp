var x{1..5} integer >= 0 <= 4;

minimize o: -4*x[1] - 4*x[2] - 2*x[3] + 0*x[4] + x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: atleast 2 {i in 1..5} (x[i] >= 3);
s.t. c3: x[1] == 4;
