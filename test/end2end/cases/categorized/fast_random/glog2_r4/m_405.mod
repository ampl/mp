var x{1..6} integer >= 0 <= 2;

minimize o: -5*x[1] - 2*x[2] + 0*x[3] - 3*x[4] - x[5] - 2*x[6];

s.t. c1: (if x[3] != 0 then x[2] else x[2]) <= 0;
s.t. c2: not ((x[2] <= 2) and (x[1] < 0));
s.t. c3: (x[2] < 1) <==> (x[5] + x[4] < 0);
s.t. c4: (x[5] != 2) <==> (x[4] + x[1] < 0);
s.t. c5: max(x[4], x[1]) - min(x[1], x[4]) == 1;
s.t. c6: not ((x[5] > 2) and (x[5] >= 0));
