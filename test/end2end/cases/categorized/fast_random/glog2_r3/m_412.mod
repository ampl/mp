var x{1..4} integer >= -2 <= 0;

maximize o: -x[1] + 0*x[2] + 2*x[3] + 0*x[4];

s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: (x[4] <= -1) or (x[4] > -2) or (x[2] < -1);
