var x{1..6} integer >= 0 <= 4;

minimize o: -x[1] - 5*x[2] - 4*x[3] + 5*x[4] - 3*x[5] + x[6];

s.t. c1: not ((x[4] == -1) and (x[2] == -1));
s.t. c2: (x[3] < 5) or (x[4] != 5) or (x[5] == -1);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: exists {i in 1..6} x[i] == 0;
