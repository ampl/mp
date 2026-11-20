var x{1..3} integer >= 0 <= 2;

maximize o: -3000*x[1] - 1000*x[2] - 2000*x[3];

s.t. c1: (x[3] > 0) and (x[3] > 0) ==> (x[2] >= 0);
s.t. c2: abs(x[1] - x[3]) <= 0;
s.t. c3: (if x[3] <= 1 then x[1] else x[3]) == 0;
s.t. c4: (x[1] >= 2) <==> (x[1] + x[2] != 2);
