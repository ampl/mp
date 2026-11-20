var x{1..5} integer >= 0 <= 2;

minimize o: 6*x[1] + x[2] + 6*x[3] - 5*x[4] - 5*x[5];

s.t. c1: (x[3] != 0) and ((x[3] <= x[2]) or not (x[3] != 0));
s.t. c2:
   (if x[3] != 0 then x[4] else 0) + (if x[3] <= x[2] then x[1] else 0) < 0;
s.t. c3: (x[3] <= x[2]) and ((x[3] <= x[2]) or not (x[3] != 0));
