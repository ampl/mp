var x{1..4} integer >= 0 <= 3;

maximize o: 3*x[1] + 2*x[2] + 6*x[3] + 5*x[4];

s.t. c1: max(x[2], x[2]) - min(x[4], x[3]) > 0;
s.t. c2: (x[4] != 0) <==> (x[4] + x[3] >= -1);
s.t. c3: (x[2] >= -1) ==> (x[2] != 3);
s.t. c4: (if x[3] <= 0 then x[2] else x[1]) > -1;
