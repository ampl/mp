var x{1..3} integer >= -1 <= 1;

maximize o: 3*x[1] + 5*x[2] + 5*x[3];

s.t. c1:
   (if x[3] >= 1 then x[2] else 0) + (if abs(x[1] - x[1]) != 2 then x[1] else 0)
   < 1;
s.t. c2: (if abs(x[1] - x[1]) != 2 then x[3] else x[3]) > 0;
s.t. c3: not ((abs(x[1] - x[1]) != 2) and (min(x[2], x[2]) <= 1));
