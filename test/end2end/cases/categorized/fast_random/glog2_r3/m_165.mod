var x{1..4} integer >= 0 <= 4;

minimize o: x[1] + 5*x[2] + 0*x[3] - 5*x[4];

s.t. c1: alldiff(x[1], x[2], x[3]);
s.t. c2: not ((x[4] <= 0) and (x[1] > 2));
s.t. c3: abs(x[2] - x[3]) == 2;
s.t. c4: exists {i in 1..4} x[i] == 4;
s.t. c5: abs(x[4] - x[3]) != 3;
s.t. c6: forall {i in 1..3} x[i] != x[i+1];
