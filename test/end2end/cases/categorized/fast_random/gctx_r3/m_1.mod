var x{1..4} integer >= -1 <= 2;

maximize o: -x[1] + x[2] - 2*x[3] - 2*x[4];

s.t. c1:
   (if abs(x[3] - x[2]) > 2 then 1 else 0)
   + (if max(x[1], x[1]) <= 0 then 1 else 0)
   + (if max(x[1], x[1]) <= 0 then 1 else 0)
   >= 3;
s.t. c2: (if max(x[4], x[4]) == 0 then x[4] else x[3]) <= 2;
s.t. c3: not (max(x[1], x[1]) <= 0);
