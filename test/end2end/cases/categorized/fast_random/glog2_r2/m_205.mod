var x{1..6} integer >= 0 <= 4;

minimize o: 4*x[1] - 5*x[2] + 6*x[3] - 5*x[4] + 4*x[5] + 2*x[6];

s.t. c1: (x[1] > 1) ==> (x[5] == 1) else (x[6] >= 4);
s.t. c2: (if x[6] == 0 then x[1] else x[6]) < 3;
s.t. c3: alldiff(x[2], x[3]);
s.t. c4: abs(x[2] - x[2]) < 3;
s.t. c5: max(x[6], x[4]) - min(x[2], x[4]) > 4;
