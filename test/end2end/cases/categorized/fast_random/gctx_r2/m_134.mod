var x{1..4} integer >= 0 <= 3;

maximize o: 5*x[1] - 3*x[2] - 3*x[3] + 4*x[4];

s.t. c1:
   (if x[4] < x[3] then 1 else 0) + (if x[1] == 0 then 1 else 0)
   + (if x[4] < x[3] then 1 else 0)
   >= 3;
s.t. c2: (if x[4] < x[3] then x[3] else x[2]) > 3;
