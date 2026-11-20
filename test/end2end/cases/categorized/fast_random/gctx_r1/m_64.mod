var x{1..5} integer >= 0 <= 2;

maximize o:
   x[1] + 2*x[2] + 3*x[3] + 5*x[4] + 2*x[5] - 3*(if x[1] > 2 then 1 else 0);

s.t. c1: not (min(x[1], x[2]) == 1);
s.t. c2: (min(x[1], x[2]) == 1) ==> (x[1] > 2);
s.t. c3: (min(x[1], x[2]) == 1) <==> (min(x[1], x[2]) == 1);
s.t. c4: (x[1] > 2) or (x[1] > 2);
