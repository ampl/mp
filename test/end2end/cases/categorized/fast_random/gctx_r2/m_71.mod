var x{1..4} integer >= 0 <= 2;

maximize o:
   -5*x[1] + 2*x[2] - 4*x[3] + 0*x[4]
   + 0*(if max(x[3], x[2]) > 0 then 1 else 0);

s.t. c1:
   (x[3] + x[2] >= 1)
   and ((abs(x[1] - x[4]) <= 2) or not (max(x[3], x[2]) > 0));
s.t. c2: (abs(x[1] - x[4]) <= 2) ==> (abs(x[1] - x[4]) <= 2);
s.t. c3: not ((x[3] + x[2] >= 1) and (x[3] + x[2] >= 1));
