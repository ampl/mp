var x{1..6} integer >= -2 <= 2;

maximize o: -4*x[1] + x[2] - x[3] - 3*x[4] - 2*x[5] + 6*x[6];

s.t. c1: ((x[6] == 0) ==> (x[3] >= 0)) and ((x[3] >= 0) ==> (x[4] >= 0));
s.t. c2: alldiff(x[1], x[6]);
s.t. c3: exists {i in 1..6} x[i] == 1;
