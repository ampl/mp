var x{1..5} integer >= -2 <= 1;

minimize o: x[1] - 4*x[2] - 2*x[3] + 0*x[4] + 6*x[5];

s.t. c1: (x[1] <= 2) and (x[5] > 2) ==> (x[2] >= 2);
s.t. c2: (if x[3] != -3 then x[2] else x[1]) > 1;
s.t. c3: (if x[2] <= -2 then x[1] else x[4]) >= -3;
