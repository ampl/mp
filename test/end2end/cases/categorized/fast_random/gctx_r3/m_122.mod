var x{1..4} integer >= -1 <= 1;

maximize o: -4*x[1] + 6*x[2] - 4*x[3] - 4*x[4];

s.t. c1: (x[4] > x[4]) and ((x[4] > x[4]) or not (x[4] > x[4]));
s.t. c2:
   (max(x[4], x[4]) <= 1)
   and ((x[3] + x[4] <= 0) or not (max(x[4], x[4]) <= 1));
