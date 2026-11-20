var x{1..5} integer >= 0 <= 2;

minimize o: -5*x[1] - 5*x[2] - 4*x[3] + 2*x[4] + 5*x[5];

s.t. c1:
   (if min(x[1], x[4]) >= 2 then 1 else 0) + (if x[4] + x[4] > 0 then 1 else 0)
   + (if x[4] + x[4] > 0 then 1 else 0)
   >= 3;
s.t. c2: (abs(x[5] - x[4]) <= 2) <==> (min(x[1], x[4]) >= 2);
s.t. c3: not ((min(x[1], x[4]) >= 2) and (min(x[1], x[4]) >= 2));
