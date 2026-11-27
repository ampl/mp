var x{1..5} integer >= 0 <= 2;

minimize o: -2*x[1] + 5*x[2] + x[3] + x[4] + 5*x[5];

s.t. c1: not ((x[5] > 2) and (x[5] <= 2));
s.t. c2: alldiff(x[5], x[1]);
s.t. c3: ((x[4] == 1) ==> (x[1] <= 2)) and ((x[1] <= 2) ==> (x[4] >= 1));
