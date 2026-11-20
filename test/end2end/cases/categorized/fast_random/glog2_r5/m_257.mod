var x{1..4} integer >= 0 <= 3;

maximize o: 17*x[1] - 999*x[2] + 1000*x[3] - 1000*x[4];

s.t. c1: ((if x[1] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[2] > 1) <==> (x[3] + x[4] < 3);
s.t. c3: (x[4] >= 1) or (x[1] == 2) or (x[1] < 1);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
