var x{1..6} integer >= 0 <= 4;

maximize o: -2*x[1] + 0*x[2] + x[3] - 5*x[4] + 0*x[5] + 0*x[6];

s.t. c1: (x[4] > 4) and (x[3] >= 5) ==> (x[1] >= 4);
s.t. c2: exists {i in 1..6} x[i] == -1;
s.t. c3: (if x[6] != 5 then x[1] else x[5]) >= 5;
