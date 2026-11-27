var x{1..6} integer >= 0 <= 3;

maximize o: 3*x[1] - 3*x[2] + 4*x[3] + x[4] + 0*x[5] - 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: abs(x[6] - x[5]) == 2;
s.t. c3: (if x[3] != 1 then x[6] else x[4]) <= 0;
s.t. c4: (x[1] != 3) <==> (x[5] + x[1] <= 1);
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
s.t. c6: ((if x[4] >= 2 then 1 else 0) + (if x[6] >= 2 then 1 else 0)) >= 1;
s.t. c7: (x[4] >= 2) <==> (x[4] + x[1] > 2);
