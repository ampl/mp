var x{1..6} integer >= -2 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (x[4] <= 1) and (x[6] <= -1) ==> (x[6] >= 1);
s.t. c3: alldiff(x[4], x[5], x[3]);
s.t. c4: ((x[6] > 1) ==> (x[6] >= 1)) and ((x[6] >= 1) ==> (x[6] >= 1));
