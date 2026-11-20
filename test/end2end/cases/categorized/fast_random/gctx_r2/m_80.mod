var x{1..3} integer >= -1 <= 1;

minimize o: 2*x[1] + x[2] - x[3] + 5*(if abs(x[1] - x[3]) <= 2 then 1 else 0);

s.t. c1: not ((min(x[3], x[1]) != 0) and (x[3] <= -1));
s.t. c2: (x[3] <= -1) <==> (abs(x[1] - x[3]) <= 2);
s.t. c3:
   (abs(x[1] - x[3]) <= 2)
   and ((abs(x[1] - x[3]) <= 2) or not (abs(x[1] - x[3]) <= 2));
