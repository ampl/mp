var x{1..3} integer >= 0 <= 3;

maximize o: -x[1] + 4*x[2] - 3*x[3];

s.t. c1: (x[3] <= 2) ==> (x[2] <= 2);
s.t. c2: (x[3] >= 2) <==> (x[2] + x[1] > 0);
s.t. c3: max(x[2], x[3]) - min(x[1], x[2]) <= 3;
s.t. c4: alldiff(x[3], x[2]);
s.t. c5: (if x[3] != 0 then x[2] else x[1]) != 3;
s.t. c6: (x[1] >= 1) and (x[2] <= 3) ==> (x[2] >= 1);
s.t. c7: forall {i in 1..2} x[i] != x[i+1];
