var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 2*x[2] - 4*x[3] + 5*x[4];

s.t. c1: atmost 1 {i in 1..4} (x[i] >= 4);
s.t. c2: alldiff(x[3], x[1], x[2]);
s.t. c3: max(x[1], x[4]) - min(x[1], x[4]) > 1;
s.t. c4: max(x[3], x[2]) - min(x[2], x[3]) == 2;
s.t. c5: exists {i in 1..4} x[i] == 0;
