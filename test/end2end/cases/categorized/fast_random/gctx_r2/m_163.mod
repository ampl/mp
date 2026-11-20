var x{1..5} integer >= 0 <= 3;

minimize o: -4*x[1] + 2*x[2] - 2*x[3] - 4*x[4] + 0*x[5];

s.t. c1:
   (if x[1] + x[5] <= 0 then x[2] else 0)
   + (if x[1] + x[5] <= 0 then x[2] else 0)
   < 0;
s.t. c2:
   (if x[2] + x[3] == 3 then 1 else 0) + (if x[1] + x[5] <= 0 then 1 else 0)
   + (if x[2] + x[3] == 3 then 1 else 0)
   >= 2;
