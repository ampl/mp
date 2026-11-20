var x{1..5} integer >= 0 <= 3;

maximize o:
   5*x[1] + x[2] - 5*x[3] + 5*x[4] - 4*x[5] - 6*(if x[3] <= 3 then 1 else 0);

s.t. c1: (x[3] <= 3) and ((x[1] <= 3) or not (x[1] <= 3));
s.t. c2:
   (if x[3] <= 3 then 1 else 0) + (if x[1] <= 3 then 1 else 0)
   + (if x[1] <= 3 then 1 else 0)
   >= 2;
