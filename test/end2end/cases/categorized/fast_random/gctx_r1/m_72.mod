var x{1..4} integer >= 0 <= 2;

minimize o: -x[1] - 3*x[2] + 6*x[3] + x[4];

s.t. c1: not ((x[2] + x[1] >= 2) and (x[2] + x[1] >= 2));
s.t. c2: not (min(x[2], x[2]) <= 1);
s.t. c3: not ((min(x[2], x[2]) <= 1) and (x[2] + x[1] >= 2));
s.t. c4:
   (min(x[2], x[2]) <= 1) ==> ((x[2] + x[1] >= 2) or (min(x[2], x[2]) <= 1));
