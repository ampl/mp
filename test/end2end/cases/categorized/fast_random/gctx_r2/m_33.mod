var x{1..5} integer >= 0 <= 2;

minimize o: 4*x[1] - 3*x[2] - 4*x[3] + 3*x[4] - 4*x[5];

s.t. c1: (max(x[5], x[4]) <= 1) or (x[5] <= 1);
s.t. c2: (max(x[5], x[5]) >= 0) <==> (max(x[5], x[5]) >= 0);
s.t. c3:
   (if x[5] <= 1 then x[2] else 0) + (if max(x[5], x[4]) <= 1 then x[1] else 0)
   > 2;
