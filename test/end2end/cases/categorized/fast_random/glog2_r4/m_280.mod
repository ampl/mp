var x{1..3} integer >= 0 <= 3;

minimize o: -3*x[1] + 3*x[2] + 2*x[3];

s.t. c1: not ((x[3] > 2) and (x[3] >= 2));
s.t. c2: ((x[2] >= 1) ==> (x[2] <= 0)) and ((x[2] <= 0) ==> (x[3] >= 1));
s.t. c3: alldiff(x[2], x[1], x[3]);
