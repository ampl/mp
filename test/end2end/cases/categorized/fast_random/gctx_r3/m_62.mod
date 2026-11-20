var x{1..5} integer >= 0 <= 3;

minimize o: -4*x[1] + 2*x[2] - 2*x[3] - 2*x[4] + 2*x[5];

s.t. c1:
   (if abs(x[1] - x[2]) <= 1 then 1 else 0) + (if x[4] >= 3 then 1 else 0)
   + (if x[4] >= 3 then 1 else 0)
   >= 2;
s.t. c2: (abs(x[1] - x[2]) <= 1) ==> ((x[4] >= 3) or (abs(x[1] - x[2]) <= 1));
