var x{1..4} integer >= 0 <= 3;

minimize o: x[1] - 3*x[2] - 2*x[3] + 2*x[4];

s.t. c1: alldiff(x[2], x[3], x[4], x[1]);
s.t. c2: abs(x[1] - x[4]) > 3;
s.t. c3: (x[3] <= 3) <==> (x[1] + x[1] >= 1);
s.t. c4: (if x[4] > 2 then x[1] else x[4]) <= 3;
s.t. c5: max(x[1], x[1]) - min(x[4], x[2]) == 0;
s.t. c6: exists {i in 1..4} x[i] == 2;
