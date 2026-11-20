var x{1..4} integer >= 0 <= 2;

minimize o: 6*x[1] - x[2] - 4*x[3] - 3*x[4];

s.t. c1: (x[3] != x[1]) ==> (abs(x[4] - x[4]) > 2) else (abs(x[4] - x[4]) > 2);
s.t. c2: (abs(x[4] - x[4]) > 2) <==> (abs(x[4] - x[4]) > 2);
s.t. c3:
   (max(x[2], x[1]) != 0)
   and ((max(x[2], x[1]) != 0) or not (max(x[2], x[1]) != 0));
