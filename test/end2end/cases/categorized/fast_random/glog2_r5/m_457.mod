var x{1..6} integer >= 0 <= 2;

minimize o: 3000*x[1] + 3000*x[2] + 3000*x[3] + 34*x[4] + 3000*x[5] - 3000*x[6];

s.t. c1: not ((x[6] != 2) and (x[4] == 0));
s.t. c2: (if x[4] != 1 then x[6] else x[1]) <= 0;
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: (x[1] == 2) and (x[2] != 0) ==> (x[3] >= 2);
