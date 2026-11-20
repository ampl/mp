var x{1..5} integer >= -2 <= 2;

maximize o: 2*x[1] + 2*x[2] + 6*x[3] + x[4] - x[5];

s.t. c1: exists {i in 1..5} x[i] == -2;
s.t. c2: alldiff(x[3], x[2], x[1]);
s.t. c3: (x[3] != -2) ==> (x[1] != -2);
s.t. c4: (x[5] == 2) ==> (x[3] != 1) else (x[2] > 0);
s.t. c5: x[3] == -2;
s.t. c6: x[5] == -2;
