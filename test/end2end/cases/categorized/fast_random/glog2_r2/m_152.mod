var x{1..5} integer >= 0 <= 4;

minimize o: 3*x[1] + 5*x[2] + 3*x[3] + 5*x[4] + 3*x[5];

s.t. c1: alldiff(x[5], x[4], x[2], x[3]);
s.t. c2: ((x[1] != 0) ==> (x[1] != 2)) and ((x[1] != 2) ==> (x[1] >= 0));
