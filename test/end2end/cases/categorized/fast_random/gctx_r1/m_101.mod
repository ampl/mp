var x{1..5} integer >= 0 <= 3;

maximize o:
   5*x[1] + 4*x[2] + 6*x[3] - 3*x[4] + x[5]
   + (if x[2] + x[4] != 3 then 1 else 0);

s.t. c1: (x[4] >= 2) or (x[4] >= 2);
s.t. c2:
   (if x[2] + x[4] != 3 then 1 else 0) + (if x[2] + x[4] != 3 then 1 else 0)
   + (if x[2] + x[4] != 3 then 1 else 0)
   >= 1;
s.t. c3: (x[2] + x[4] != 3) and ((x[4] >= 2) or not (x[4] >= 2));
s.t. c4: (x[4] >= 2) or (x[4] >= 2);
