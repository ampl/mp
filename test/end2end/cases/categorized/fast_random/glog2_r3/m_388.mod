var x{1..7} integer >= 0 <= 1;

maximize o: 5*x[1] - 5*x[2] - 3*x[3] + 2*x[4] + 6*x[5] + 3*x[6] + 3*x[7];

s.t. c1: (x[7] >= 1) and (x[4] < 0) ==> (x[6] >= 1);
s.t. c2: (x[6] > 1) or (x[5] > 0) or (x[4] == 0);
s.t. c3: alldiff(x[7], x[2], x[3]);
