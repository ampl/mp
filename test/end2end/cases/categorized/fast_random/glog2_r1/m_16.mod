var x{1..6} integer >= -2 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: abs(x[6] - x[5]) == 3;
s.t. c2: (if x[1] <= 0 then x[4] else x[4]) == 1;
s.t. c3: (x[2] >= 0) or (x[5] > 2) or (x[4] < -1);
