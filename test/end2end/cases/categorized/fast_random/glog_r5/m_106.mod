var x{1..5} integer >= 0 <= 3;

minimize o: -x[1] - 2*x[2] - 4*x[3] + 3*x[4] - 4*x[5];

s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: alldiff(x[5], x[2]);
s.t. c3: (x[5] > 2) ==> (x[4] >= 3);
