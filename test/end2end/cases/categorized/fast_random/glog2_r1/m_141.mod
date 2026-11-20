var x{1..3} integer >= 0 <= 3;

maximize o: 5*x[1] - 2*x[2] + 4*x[3];

s.t. c1: abs(x[3] - x[2]) < 0;
s.t. c2: exists {i in 1..3} x[i] == -1;
