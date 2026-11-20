var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] + 6*x[2] - 5*x[3];

s.t. c1: abs(x[2] - x[1]) > 2;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (if x[2] > 4 then x[3] else x[3]) < 4;
s.t. c4: alldiff(x[1], x[2], x[3]);
