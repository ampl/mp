var x{1..6} integer >= 0 <= 3;

maximize o: -5*x[1] + 6*x[2] - 2*x[3] - 2*x[4] + 6*x[5] - 3*x[6];

s.t. c1: (if x[3] == 1 then x[3] else x[6]) >= 2;
s.t. c2: abs(x[6] - x[1]) >= 3;
s.t. c3: (x[1] >= 1) <==> (x[3] + x[5] <= 2);
s.t. c4: not ((x[3] > 0) and (x[2] > 0));
s.t. c5: exists {i in 1..6} x[i] == 3;
s.t. c6: (x[5] > 1) or (x[6] == 0) or (x[4] >= 0);
s.t. c7: (x[6] != 1) ==> (x[5] > 1) else (x[5] == 1);
