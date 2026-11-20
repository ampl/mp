var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] + 2*x[2] - 5*x[3] + 5*x[4];

s.t. c1: (if x[2] > 3 then x[2] else x[4]) == 0;
s.t. c2: (x[2] != 0) and (x[4] <= 2) ==> (x[2] >= 0);
