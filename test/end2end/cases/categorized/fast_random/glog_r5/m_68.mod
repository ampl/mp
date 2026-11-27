var x{1..4} integer >= 0 <= 3;

maximize o: 2*x[1] + 0*x[2] - 5*x[3] - 5*x[4];

s.t. c1: (x[4] > 0) or (x[4] <= 1) or (x[4] >= 2);
s.t. c2: (if x[1] < 2 then x[3] else x[2]) > 2;
s.t. c3: (x[3] <= 1) ==> (x[4] >= 2) else (x[4] < 2);
s.t. c4: (x[2] < 3) and (x[2] > 0) ==> (x[4] >= 3);
