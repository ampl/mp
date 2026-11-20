var x{1..4} integer >= 0 <= 2;

maximize o:
   -2*x[1] - x[2] - x[3] + x[4] - 3*(if abs(x[1] - x[3]) == 2 then 1 else 0);

s.t. c1: not (abs(x[1] - x[3]) == 2);
s.t. c2: (abs(x[1] - x[3]) == 2) ==> (x[3] <= x[1]);
s.t. c3: (abs(x[1] - x[3]) == 2) <==> (x[3] <= x[1]);
