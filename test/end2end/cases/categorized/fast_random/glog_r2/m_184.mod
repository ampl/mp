var x{1..6} integer >= 0 <= 2;

minimize o: 4*x[1] + x[2] + 5*x[3] + 5*x[4] + 2*x[5] + 5*x[6];

s.t. c1: (if x[6] >= 1 then x[4] else x[5]) != 1;
s.t. c2: exists {i in 1..6} x[i] == 1;
