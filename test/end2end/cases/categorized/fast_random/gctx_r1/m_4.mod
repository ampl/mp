var x{1..3} integer >= -1 <= 1;

minimize o: -2*x[1] - 5*x[2] - 4*x[3];

s.t. c1: (x[1] + x[3] < 0) <==> (x[1] + x[3] < 0);
s.t. c2:
   (if abs(x[2] - x[2]) < 2 then x[2] else 0)
   + (if abs(x[2] - x[2]) < 2 then x[1] else 0)
   >= 0;
s.t. c3: (x[1] + x[3] < 0) ==> (x[1] + x[3] < 0) else (abs(x[2] - x[2]) < 2);
