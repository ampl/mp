var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] - 4*x[2] - 2*x[3] + 6*x[4];

s.t. c1:
   (min(x[4], x[1]) != 3)
   and ((min(x[4], x[1]) != 3) or not (abs(x[3] - x[1]) != 2));
s.t. c2: (min(x[4], x[1]) != 3) or (max(x[2], x[1]) > 3);
