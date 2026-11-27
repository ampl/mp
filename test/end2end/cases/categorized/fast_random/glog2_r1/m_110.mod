var x{1..5} integer >= 0 <= 4;

minimize o: 3*x[1] - 4*x[2] + 2*x[3] + 4*x[4] - 2*x[5];

s.t. c1: (if x[3] != 1 then x[1] else x[5]) < 4;
s.t. c2: (x[4] > 2) ==> (x[2] < 2);
s.t. c3: not ((x[3] < 0) and (x[4] >= 3));
