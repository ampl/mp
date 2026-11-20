var x{1..7} integer >= 0 <= 1;

maximize o: 0*x[1] + 6*x[2] - 2*x[3] + 2*x[4] - x[5] - 4*x[6] + 0*x[7];

s.t. c1: max(x[3], x[5]) - min(x[5], x[2]) == 0;
s.t. c2: exists {i in 1..7} x[i] == 1;
s.t. c3: alldiff(x[4], x[3], x[7], x[5]);
s.t. c4: (x[1] != 1) ==> (x[4] != 0);
