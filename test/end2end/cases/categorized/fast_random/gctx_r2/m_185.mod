var x{1..3} integer >= -1 <= 1;

minimize o: 5*x[1] + 0*x[2] + 3*x[3] - 2*(if x[2] <= -1 then 1 else 0);

s.t. c1:
   (x[2] <= -1) and ((min(x[2], x[2]) <= 0) or not (min(x[2], x[2]) <= 0));
s.t. c2: (if x[2] <= -1 then x[3] else x[3]) > 0;
