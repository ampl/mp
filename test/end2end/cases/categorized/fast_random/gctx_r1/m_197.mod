var x{1..5} integer >= -1 <= 2;

maximize o: 5*x[1] + 4*x[2] + 5*x[3] + 4*x[4] - x[5];

s.t. c1:
   (if x[5] + x[1] < 0 then 1 else 0) + (if x[5] + x[1] < 0 then 1 else 0)
   + (if x[2] != x[4] then 1 else 0)
   >= 2;
s.t. c2: (x[5] + x[1] < 0) and ((x[5] + x[1] < 0) or not (x[4] != x[3]));
