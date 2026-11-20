var x{1..3} integer >= -1 <= 2;

maximize o: -4*x[1] + 2*x[2] + 3*x[3];

s.t. c1:
   (if abs(x[1] - x[1]) < 1 then 1 else 0)
   + (if abs(x[1] - x[1]) < 1 then 1 else 0)
   + (if x[3] + x[1] >= 2 then 1 else 0)
   >= 3;
s.t. c2:
   (abs(x[1] - x[1]) < 1)
   ==> ((min(x[3], x[3]) == 2) or (abs(x[1] - x[1]) < 1));
s.t. c3:
   (abs(x[1] - x[1]) < 1)
   ==> (min(x[3], x[3]) == 2)
   else (abs(x[1] - x[1]) < 1);
s.t. c4:
   (abs(x[1] - x[1]) < 1)
   and ((x[3] + x[1] >= 2) or not (abs(x[1] - x[1]) < 1));
