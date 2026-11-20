var x{1..5} integer >= -2 <= 0;

minimize o: 5*x[1] + 5*x[2] + 4*x[3] + 6*x[4] + 3*x[5];

s.t. c1: (x[3] < -1) ==> (x[1] < -2);
s.t. c2: not ((x[2] >= -1) and (x[5] < -2));
s.t. c3: exists {i in 1..5} x[i] == -1;
