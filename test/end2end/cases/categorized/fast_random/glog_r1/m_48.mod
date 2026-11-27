var x{1..6} integer >= 0 <= 4;

maximize o: -x[1] + 4*x[2] + 3*x[3] + 2*x[4] + x[5] - 4*x[6];

s.t. c1: abs(x[1] - x[2]) == 2;
s.t. c2: (if x[5] != 0 then x[5] else x[5]) < 1;
s.t. c3: not ((x[6] < 3) and (x[3] == 4));
s.t. c4: exists {i in 1..6} x[i] == 3;
