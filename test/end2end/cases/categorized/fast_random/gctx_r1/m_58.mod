var x{1..5} integer >= -1 <= 1;

maximize o:
   2*x[1] - x[2] + 5*x[3] + 0*x[4] + 5*x[5]
   + 0*(if abs(x[3] - x[2]) == 1 then 1 else 0);

s.t. c1: not ((max(x[1], x[4]) > -1) and (max(x[1], x[4]) > -1));
s.t. c2: (abs(x[3] - x[2]) == 1) ==> (max(x[3], x[2]) != -1);
