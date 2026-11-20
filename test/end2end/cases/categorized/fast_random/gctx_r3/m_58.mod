var x{1..3} integer >= 0 <= 3;

maximize o: 0*x[1] - 2*x[2] + 6*x[3] + 3*(if x[1] + x[1] > 1 then 1 else 0);

s.t. c1: (min(x[1], x[3]) < 3) <==> (x[1] + x[1] > 1);
s.t. c2: (x[2] == 2) ==> (x[2] == 2) else (x[2] == 2);
