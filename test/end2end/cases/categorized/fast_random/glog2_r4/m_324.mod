var x{1..5} integer >= -2 <= 2;

maximize o: 2*x[1] - 4*x[2] - x[3] - 5*x[4] - 3*x[5];

s.t. c1: not ((x[2] == 2) and (x[2] >= 2));
s.t. c2: (x[5] <= 2) ==> (x[3] <= -1) else (x[5] != -2);
s.t. c3: (x[1] > 1) and (x[5] <= -1) ==> (x[1] >= 1);
s.t. c4: alldiff(x[5], x[3], x[4], x[2]);
