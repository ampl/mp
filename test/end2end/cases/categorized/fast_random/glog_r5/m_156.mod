var x{1..6} integer >= 0 <= 3;

maximize o: 4*x[1] - x[2] - 2*x[3] - x[4] + x[5] + 6*x[6];

s.t. c1: (x[2] == 3) ==> (x[1] != 1);
s.t. c2: (x[5] == 3) <==> (x[6] + x[1] < 2);
s.t. c3: (if x[2] >= 1 then x[2] else x[1]) < 2;
