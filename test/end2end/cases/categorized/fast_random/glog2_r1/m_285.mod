var x{1..6} integer >= -2 <= 2;

maximize o: x[1] - 3*x[2] + 6*x[3] - 4*x[4] - x[5] + 2*x[6];

s.t. c1: max(x[3], x[6]) - min(x[6], x[1]) < 2;
s.t. c2: (x[1] > 0) and (x[3] >= -1) ==> (x[4] >= 0);
s.t. c3: (if x[6] != 0 then x[3] else x[6]) == 1;
s.t. c4: abs(x[2] - x[1]) >= 4;
s.t. c5: exists {i in 1..6} x[i] == -1;
s.t. c6: atleast 3 {i in 1..6} (x[i] >= 1);
