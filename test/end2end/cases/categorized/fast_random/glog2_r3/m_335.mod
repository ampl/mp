var x{1..6} integer >= 0 <= 3;

maximize o: x[1] - 4*x[2] + 4*x[3] - 2*x[4] + 3*x[5] - 4*x[6];

s.t. c1: (x[1] < 3) <==> (x[5] + x[3] == 1);
s.t. c2: alldiff(x[5], x[4], x[1], x[6]);
s.t. c3: ((if x[1] < 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c4: (x[5] < 2) and (x[5] >= 0) ==> (x[5] >= 2);
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
s.t. c6: (x[2] > 0) ==> (x[4] > 0);
