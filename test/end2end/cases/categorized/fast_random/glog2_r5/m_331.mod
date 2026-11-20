var x{1..6} integer >= -2 <= 0;

maximize o: x[1] + 5*x[2] - 4*x[3] + x[4] + 3*x[5] + 4*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: ((x[6] >= 1) ==> (x[6] <= 1)) and ((x[6] <= 1) ==> (x[5] >= 1));
s.t. c3: exists {i in 1..6} x[i] == -2;
s.t. c4: (if x[2] == 1 then x[5] else x[6]) == -3;
