var x{1..4} integer >= 0 <= 2;

minimize o: x[1] + 4*x[2] - 3*x[3] + 0*x[4];

s.t. c1: (x[4] > 2) ==> (x[1] > 2);
s.t. c2: not ((x[3] < 0) and (x[4] > 2));
s.t. c3: (x[3] == 2) <==> (x[1] + x[4] <= 1);
s.t. c4: abs(x[2] - x[1]) < 2;
s.t. c5: not ((x[4] >= 1) and (x[2] < 0));
s.t. c6: forall {i in 1..3} x[i] != x[i+1];
