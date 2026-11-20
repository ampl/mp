var x{1..6} integer >= 0 <= 4;

minimize o: -3*x[1] + 5*x[2] - 4*x[3] + 3*x[4] + 2*x[5] + 0*x[6];

s.t. c1: not ((x[4] < 3) and (x[3] >= 3));
s.t. c2: (if x[6] <= 0 then x[2] else x[3]) >= 0;
s.t. c3: (x[2] < 3) and (x[6] != 0) ==> (x[3] >= 3);
