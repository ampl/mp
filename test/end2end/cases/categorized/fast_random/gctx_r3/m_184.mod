var x{1..3} integer >= 0 <= 3;

minimize o: 3*x[1] - 4*x[2] - 3*x[3];

s.t. c1:
   (if x[1] < x[2] then 1 else 0) + (if x[1] < x[2] then 1 else 0)
   + (if max(x[1], x[1]) != 2 then 1 else 0)
   >= 3;
s.t. c2: (x[1] < x[2]) ==> ((x[1] < x[2]) or (x[1] < x[2]));
