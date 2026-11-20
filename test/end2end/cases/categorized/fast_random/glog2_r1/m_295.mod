var x{1..6} integer >= 0 <= 3;

maximize o: 5*x[1] + x[2] + 4*x[3] + 2*x[4] + 3*x[5] + 6*x[6];

s.t. c1: atmost 1 {i in 1..6} (x[i] >= 0);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: numberof 3 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 3;
s.t. c4: alldiff(x[2], x[4], x[6]);
s.t. c5: atleast 1 {i in 1..6} (x[i] >= 2);
s.t. c6: (x[1] != 2) <==> (x[3] + x[4] == 1);
s.t. c7: ((if x[5] != 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
