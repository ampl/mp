var x{1..5} integer >= 0 <= 2;

maximize o: 5*x[1] - 4*x[2] + 3*x[3] + 0*x[4] - 5*x[5];

s.t. c1: (x[3] < 0) and (x[5] != 0) ==> (x[3] >= 0);
s.t. c2: alldiff(x[1], x[5], x[2]);
