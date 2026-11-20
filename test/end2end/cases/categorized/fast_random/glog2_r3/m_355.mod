var x{1..5} integer >= -2 <= 2;

minimize o: -x[1] + 0*x[2] - 5*x[3] + 4*x[4] - x[5];

s.t. c1: (x[4] != -1) and (x[1] != -2) ==> (x[3] >= -1);
s.t. c2: abs(x[5] - x[5]) >= 3;
s.t. c3: (x[5] <= -2) ==> (x[5] > -1);
s.t. c4: ((if x[2] > -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c5: (x[3] <= 0) ==> (x[5] != -1);
s.t. c6: max(x[2], x[3]) - min(x[1], x[2]) != -2;
