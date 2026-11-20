var x{1..4} integer >= -1 <= 1;

maximize o: -5*x[1] + 5*x[2] + 2*x[3] - x[4];

s.t. c1:
   (if x[1] < 1 then 1 else 0) + (if x[1] < 1 then 1 else 0)
   + (if x[4] >= x[2] then 1 else 0)
   >= 1;
s.t. c2: (if x[4] >= x[2] then x[1] else x[1]) >= -1;
