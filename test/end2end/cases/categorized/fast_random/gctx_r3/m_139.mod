var x{1..4} integer >= 0 <= 2;

minimize o:
   -2*x[1] + 2*x[2] - 3*x[3] + 3*x[4]
   - 5*(if min(x[2], x[2]) < 1 then 1 else 0);

s.t. c1: (x[1] >= 2) <==> (x[1] >= 2);
s.t. c2: (x[1] >= 2) or (min(x[1], x[4]) >= 1);
s.t. c3:
   (min(x[1], x[4]) >= 1)
   and ((min(x[2], x[2]) < 1) or not (min(x[2], x[2]) < 1));
s.t. c4: (min(x[2], x[2]) < 1) and ((x[4] + x[4] <= 1) or not (x[1] >= 2));
