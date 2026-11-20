var x{1..3} integer >= -1 <= 1;

maximize o: -3*x[1] - 5*x[2] + 5*x[3];

s.t. c1: not ((abs(x[3] - x[3]) >= 1) and (x[2] <= x[1]));
s.t. c2:
   (if x[2] <= x[1] then x[2] else 0) + (if x[2] <= x[1] then x[2] else 0) >= 0;
s.t. c3:
   (x[2] <= x[1]) and ((abs(x[3] - x[3]) >= 1) or not (abs(x[3] - x[3]) >= 1));
