var x{1..4} integer >= 0 <= 3;

maximize o: 4*x[1] + 5*x[2] + 5*x[3] + 5*x[4];

s.t. c1: (x[4] < 0) ==> (x[3] > 2);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: ((if x[2] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[3], x[1], x[4], x[2]);
s.t. c5: ((if x[2] <= 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c6: numberof 0 in (x[1], x[2], x[3], x[4]) < 3;
