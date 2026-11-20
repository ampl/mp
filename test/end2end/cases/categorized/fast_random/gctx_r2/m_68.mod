var x{1..3} integer >= -1 <= 2;

minimize o: 0*x[1] + 6*x[2] + 6*x[3] + 0*(if x[2] + x[2] != 2 then 1 else 0);

s.t. c1:
   (if x[2] <= 2 then 1 else 0) + (if x[2] >= 0 then 1 else 0)
   + (if x[2] <= 2 then 1 else 0)
   >= 3;
s.t. c2: (x[3] >= x[3]) ==> ((x[2] + x[2] != 2) or (x[2] >= 0));
s.t. c3: (if x[2] >= 0 then x[2] else 0) + (if x[2] >= 0 then x[2] else 0) >= 2;
