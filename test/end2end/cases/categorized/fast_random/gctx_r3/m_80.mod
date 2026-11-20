var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] - 4*x[2] + 6*x[3] - x[4];

s.t. c1: (abs(x[3] - x[2]) < 0) <==> (min(x[2], x[4]) == 1);
s.t. c2: (x[1] != 0) or (x[1] > 3);
s.t. c3:
   (if x[1] > 3 then x[1] else 0) + (if abs(x[3] - x[2]) < 0 then x[4] else 0)
   < 0;
s.t. c4: (x[1] > 3) ==> (min(x[2], x[4]) == 1) else (x[1] > 3);
