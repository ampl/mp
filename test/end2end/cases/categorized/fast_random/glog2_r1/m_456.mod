var x{1..5} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: (x[4] == 0) and (x[2] < 2) ==> (x[4] >= 0);
s.t. c2: alldiff(x[3], x[2], x[4], x[1]);
s.t. c3: not ((x[5] >= 3) and (x[4] == 0));
