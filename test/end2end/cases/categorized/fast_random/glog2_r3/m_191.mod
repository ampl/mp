var x{1..6} integer >= -2 <= 2;

maximize o: -x[1] + 3*x[2] + 4*x[3] + 6*x[4] - 3*x[5] + 5*x[6];

s.t. c1: (x[2] >= -2) ==> (x[4] < 2);
s.t. c2: (if x[1] >= 3 then x[3] else x[3]) >= 2;
s.t. c3: not ((x[3] == -2) and (x[2] != 2));
s.t. c4: (x[5] != 3) or (x[1] < -2) or (x[4] >= 2);
