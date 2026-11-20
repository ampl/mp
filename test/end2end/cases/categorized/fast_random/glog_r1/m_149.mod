var x{1..6} integer >= -2 <= 2;

maximize o: -2*x[1] + 2*x[2] + 2*x[3] - x[4] - 4*x[5] - x[6];

s.t. c1: atleast 3 {i in 1..6} (x[i] >= -2);
s.t. c2: (if x[1] == -2 then x[4] else x[2]) <= -2;
