var x{1..6} integer >= 0 <= 2;

maximize o: -x[1] + 2*x[2] + 4*x[3] + 4*x[4] + 5*x[5] - 2*x[6];

s.t. c1: max(x[5], x[4]) - min(x[6], x[2]) != 1;
s.t. c2: (x[4] == 2) ==> (x[1] <= 0);
s.t. c3: max(x[6], x[5]) - min(x[5], x[3]) > 1;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
