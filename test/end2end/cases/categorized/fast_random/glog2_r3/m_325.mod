var x{1..6} integer >= 0 <= 2;

maximize o: -3*x[1] - x[2] + 4*x[3] - 4*x[4] + 0*x[5] + 6*x[6];

s.t. c1: (x[3] > 0) <==> (x[1] + x[1] <= 1);
s.t. c2: (x[1] > 2) and (x[3] == 2) ==> (x[4] >= 2);
s.t. c3: (x[2] <= 1) <==> (x[3] + x[6] > 0);
s.t. c4: max(x[5], x[4]) - min(x[3], x[1]) >= 1;
s.t. c5: not ((x[1] <= 2) and (x[3] != 1));
s.t. c6: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 0;
s.t. c7: forall {i in 1..5} x[i] != x[i+1];
