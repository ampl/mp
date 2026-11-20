var x{1..5} integer >= 0 <= 2;

minimize o:
   -x[1] + 0*x[2] + x[3] - 2*x[4] - 5*x[5]
   - 6*(if x[3] + x[5] > 1 then 1 else 0);

s.t. c1: (if x[3] + x[5] > 1 then x[5] else x[1]) < 2;
s.t. c2: not (x[5] != 0);
s.t. c3: (if max(x[5], x[4]) == 2 then x[2] else x[5]) <= 1;
