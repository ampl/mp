var x{1..3} integer >= 0 <= 3;

minimize o: -2*x[1] + 4*x[2] - 5*x[3];

s.t. c1: not ((x[1] > 0) and (x[1] == 3));
s.t. c2: (x[1] < 2) and (x[2] != 0) ==> (x[1] >= 2);
s.t. c3: alldiff(x[2], x[3]);
