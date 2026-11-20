var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] - 2*x[2] - x[3] + x[4];

s.t. c1:
   (if abs(x[2] - x[1]) != 1 then 1 else 0)
   + (if abs(x[2] - x[1]) != 1 then 1 else 0)
   + (if abs(x[3] - x[4]) != 0 then 1 else 0)
   >= 2;
s.t. c2:
   (abs(x[2] - x[1]) != 1) ==> (abs(x[2] - x[1]) != 1) else (x[3] + x[4] < 2);
