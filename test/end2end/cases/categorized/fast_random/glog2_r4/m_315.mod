var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] - 3*x[2] + 6*x[3] + 6*x[4];

s.t. c1: max(x[2], x[3]) - min(x[3], x[4]) <= 0;
s.t. c2: (x[2] == 3) and (x[3] <= 0) ==> (x[3] >= 3);
s.t. c3: (x[3] >= 2) ==> (x[4] < 3);
s.t. c4: ((if x[3] <= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c5: not ((x[3] < 3) and (x[3] < 1));
s.t. c6: ((x[1] == 1) ==> (x[1] < 3)) and ((x[1] < 3) ==> (x[3] >= 1));
