var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] + 3*x[2] + 6*x[3] - (if x[3] < x[2] then 1 else 0);

s.t. c1: not ((x[3] < x[2]) and (x[1] < 1));
s.t. c2: not (x[3] < x[2]);
s.t. c3:
   (if x[3] < x[2] then x[1] else 0) + (if x[3] < x[2] then x[2] else 0) < 2;
