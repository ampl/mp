var x{1..3} integer >= -1 <= 2;

minimize o: 5*x[1] + x[2] + x[3];

s.t. c1: (x[2] >= 1) <==> (x[2] >= 1);
s.t. c2: not (x[2] >= 1);
s.t. c3:
   (if abs(x[3] - x[1]) == 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)
   + (if abs(x[3] - x[1]) == 1 then 1 else 0)
   >= 1;
