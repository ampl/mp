var x{1..3} integer >= 0 <= 2;

maximize o: -3*x[1] + x[2] + 5*x[3];

s.t. c1: (x[1] != 0) and (x[1] < 2) ==> (x[1] >= 0);
s.t. c2: alldiff(x[3], x[1], x[2]);
s.t. c3: ((x[3] < 1) ==> (x[3] >= 2)) and ((x[3] >= 2) ==> (x[3] >= 1));
