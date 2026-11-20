var x{1..3} integer >= -2 <= 2;

minimize o: 5*x[1] + 6*x[2] + 4*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[1], x[3]);
s.t. c3: max(x[2], x[2]) - min(x[3], x[3]) == 0;
s.t. c4: exists {i in 1..3} x[i] == -2;
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: exists {i in 1..3} x[i] == -2;
