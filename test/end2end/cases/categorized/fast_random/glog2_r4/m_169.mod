var x{1..6} integer >= -2 <= 2;
minimize o: (5)*x[1] + (6)*x[2] + (1)*x[3] + (-4)*x[4] + (1)*x[5] + (4)*x[6];
s.t. c1: max(x[1], x[3]) - min(x[6], x[1]) <= -1;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: max(x[3], x[1]) - min(x[2], x[6]) < 2;
s.t. c4: max(x[1], x[3]) - min(x[6], x[1]) <= -1;
s.t. c5: max(x[3], x[1]) - min(x[2], x[6]) < 2;
