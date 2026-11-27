var x{1..4} integer >= 0 <= 2;

minimize o: 3*x[1] - 5*x[2] - 4*x[3] - 3*x[4];

s.t. c1: (x[4] != 2) ==> (x[4] <= 2);
s.t. c2: ((if x[1] < 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[1], x[3], x[2]);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: (if x[2] < 1 then x[2] else x[3]) != 1;
s.t. c6: atmost 1 {i in 1..4} (x[i] >= 1);
