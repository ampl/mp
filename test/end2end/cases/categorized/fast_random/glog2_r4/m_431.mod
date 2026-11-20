var x{1..5} integer >= -2 <= 1;

minimize o: x[1] - x[2] + 4*x[3] - 5*x[4] + x[5];

s.t. c1: alldiff(x[3], x[1], x[4]);
s.t. c2: exists {i in 1..5} x[i] == -3;
s.t. c3: (x[4] < 1) ==> (x[2] == -3) else (x[4] < 1);
