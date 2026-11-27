var x{1..6} integer >= -2 <= 2;

maximize o: 0*x[1] - 2*x[2] + 2*x[3] + x[4] + 3*x[5] + 4*x[6];

s.t. c1: alldiff(x[1], x[2]);
s.t. c2: (x[5] < -2) ==> (x[5] > 1) else (x[2] < 0);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: alldiff(x[3], x[2], x[5]);
s.t. c5: (if x[2] <= 1 then x[4] else x[3]) == 1;
s.t. c6: exists {i in 1..6} x[i] == -1;
s.t. c7: numberof -2 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 2;
