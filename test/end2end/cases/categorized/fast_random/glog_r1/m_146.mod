var x{1..5} integer >= 0 <= 2;

maximize o: -x[1] + 3*x[2] - 2*x[3] + x[4] - 4*x[5];

s.t. c1: abs(x[1] - x[5]) <= 2;
s.t. c2: exists {i in 1..5} x[i] == 2;
