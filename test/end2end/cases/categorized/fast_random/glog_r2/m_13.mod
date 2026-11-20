var x{1..6} integer >= -2 <= 2;

maximize o: 3*x[1] + 6*x[2] - 4*x[3] + 2*x[4] - 5*x[5] + 3*x[6];

s.t. c1: (if x[6] != 2 then x[6] else x[6]) <= 2;
s.t. c2: (x[6] > -2) <==> (x[3] + x[3] >= 2);
