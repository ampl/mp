var x{1..4} integer >= -2 <= 0;

maximize o: 0*x[1] - 3*x[2] + 4*x[3] + 4*x[4];

s.t. c1: max(x[1], x[4]) - min(x[2], x[4]) <= -2;
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: (if x[2] == 0 then x[2] else x[4]) == -1;
s.t. c4: (if x[4] == -2 then x[3] else x[1]) != -2;
s.t. c5: atleast 3 {i in 1..4} (x[i] >= -2);
s.t. c6: abs(x[2] - x[2]) > 1;
