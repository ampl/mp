var x{1..5} integer >= 0 <= 2;

minimize o: 5*x[1] + 6*x[2] + 4*x[3] - 4*x[4] + 5*x[5];

s.t. c1: max(x[1], x[1]) - min(x[4], x[4]) < 2;
s.t. c2: abs(x[5] - x[5]) < 2;
s.t. c3: exists {i in 1..5} x[i] == 2;
s.t. c4: (if x[3] <= 0 then x[4] else x[1]) < 0;
s.t. c5: max(x[1], x[1]) - min(x[4], x[4]) < 2;
s.t. c6: (if x[3] <= 0 then x[4] else x[1]) < 0;
