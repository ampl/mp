var x{1..3} integer >= 0 <= 2;

maximize o: 0*x[1] + x[2] - 2*x[3] + 0*(if abs(x[2] - x[2]) < 2 then 1 else 0);

s.t. c1: (if abs(x[2] - x[2]) < 2 then x[2] else x[3]) > 2;
s.t. c2:
   (if abs(x[2] - x[2]) < 2 then 1 else 0)
   + (if abs(x[2] - x[2]) < 2 then 1 else 0)
   + (if abs(x[2] - x[2]) < 2 then 1 else 0)
   >= 3;
s.t. c3:
   (abs(x[3] - x[1]) > 1)
   and ((abs(x[2] - x[2]) < 2) or not (max(x[2], x[1]) == 0));
