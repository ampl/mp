var x{1..4} integer >= -2 <= 1;

minimize o: -5*x[1] + 3*x[2] + 0*x[3] - 4*x[4];

s.t. c1: (if x[1] <= 1 then x[3] else x[1]) >= -2;
s.t. c2: not ((x[4] > -3) and (x[3] <= -3));
s.t. c3: exists {i in 1..4} x[i] == 2;
