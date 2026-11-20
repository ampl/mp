var x{1..5} integer >= -3 <= 2;

maximize o: -4*x[1] + 5*x[2] + 6*x[3] + 6*x[4] + x[5];

s.t. c1: exists {i in 1..5} x[i] == -1;
s.t. c2: (x[3] < 2) ==> (x[3] < -2) else (x[5] != 1);
s.t. c3: alldiff(x[1], x[2], x[5]);
