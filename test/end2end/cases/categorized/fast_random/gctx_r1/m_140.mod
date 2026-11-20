var x{1..3} integer >= 0 <= 2;

maximize o: -x[1] - 2*x[2] + 5*x[3];

s.t. c1:
   (if x[3] + x[3] >= 2 then x[2] else 0)
   + (if x[3] + x[3] >= 2 then x[2] else 0)
   <= 2;
s.t. c2: not ((min(x[3], x[3]) != 2) and (x[3] + x[3] >= 2));
s.t. c3: (min(x[3], x[3]) != 2) or (min(x[3], x[3]) != 2);
