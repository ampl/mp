var x{1..5} integer >= -2 <= 0;

minimize o: 4*x[1] + 6*x[2] + 2*x[3] + 2*x[4] - 3*x[5];

s.t. c1: (if x[1] > -2 then x[1] else x[2]) != -1;
s.t. c2: (x[5] <= -1) and (x[3] > 0) ==> (x[3] >= -1);
s.t. c3: max(x[2], x[3]) - min(x[2], x[2]) == -2;
s.t. c4: (x[2] >= -2) or (x[2] > -2) or (x[3] >= -1);
