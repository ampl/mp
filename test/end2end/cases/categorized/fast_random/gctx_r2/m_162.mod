var x{1..4} integer >= 0 <= 3;

maximize o:
   -3*x[1] + 0*x[2] - 3*x[3] - 2*x[4]
   - 2*(if min(x[1], x[3]) != 3 then 1 else 0);

s.t. c1: (x[1] >= 0) or (x[3] <= 3);
s.t. c2: (x[1] >= 0) ==> ((min(x[1], x[3]) != 3) or (min(x[1], x[3]) != 3));
s.t. c3: (min(x[1], x[3]) != 3) or (x[1] >= 0);
