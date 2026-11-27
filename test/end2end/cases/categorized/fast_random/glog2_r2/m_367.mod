var x{1..5} integer >= 0 <= 3;

maximize o: 34*x[1] - 3000*x[2] + 3000*x[3] - 3000*x[4] + 51*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (if x[3] <= 1 then x[4] else x[5]) >= 0;
s.t. c3: (x[4] < 3) and (x[2] != 1) ==> (x[5] >= 3);
s.t. c4: (x[5] >= 0) <==> (x[5] + x[5] != 0);
