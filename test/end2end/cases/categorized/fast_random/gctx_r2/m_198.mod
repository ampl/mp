var x{1..4} integer >= -1 <= 1;

minimize o: x[1] - 4*x[2] + x[3] + 6*x[4];

s.t. c1: not (x[2] + x[3] >= 1);
s.t. c2:
   (if x[1] < 0 then x[2] else 0) + (if x[2] + x[3] >= 1 then x[2] else 0) <= 1;
s.t. c3: (if x[2] + x[3] >= 1 then x[1] else x[4]) == -1;
s.t. c4: (if x[2] + x[3] >= 1 then x[1] else x[3]) > 1;
