var x{1..5} integer >= -2 <= 2;

minimize o: x[1] - 5*x[2] - 4*x[3] - 5*x[4] + 5*x[5];

s.t. c1: (x[2] == 1) and (x[3] <= -2) ==> (x[3] >= 1);
s.t. c2: ((if x[4] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: max(x[2], x[5]) - min(x[2], x[4]) == -1;
s.t. c4: (x[1] >= -1) ==> (x[2] <= 1);
