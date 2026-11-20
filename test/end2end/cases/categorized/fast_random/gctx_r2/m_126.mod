var x{1..4} integer >= -1 <= 1;

minimize o: 3*x[1] + 0*x[2] + 0*x[3] - 5*x[4];

s.t. c1: (x[2] != -1) <==> (abs(x[3] - x[3]) == 1);
s.t. c2:
   (if x[2] != -1 then x[2] else 0) + (if min(x[1], x[2]) >= 1 then x[3] else 0)
   <= -1;
s.t. c3:
   (if x[3] + x[3] == 0 then 1 else 0) + (if x[3] + x[3] == 0 then 1 else 0)
   + (if min(x[1], x[2]) >= 1 then 1 else 0)
   >= 3;
