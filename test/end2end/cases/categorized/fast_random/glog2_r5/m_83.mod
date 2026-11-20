var x{1..4} integer >= 0 <= 4;

maximize o: -3*x[1] - 2*x[2] + 5*x[3] + 3*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: abs(x[2] - x[1]) == 4;
s.t. c3: alldiff(x[3], x[2], x[4], x[1]);
s.t. c4: x[4] == 1;
s.t. c5: x[3] == 4;
