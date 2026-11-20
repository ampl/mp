var x{1..3} integer >= 0 <= 3;

minimize o: -2*x[1] + 0*x[2] + 5*x[3];

s.t. c1:
   (if abs(x[2] - x[2]) == 3 then x[1] else 0)
   + (if abs(x[2] - x[2]) == 3 then x[2] else 0)
   < 1;
s.t. c2: not ((abs(x[2] - x[2]) == 3) and (x[1] == x[3]));
s.t. c3: not ((x[1] == x[3]) and (x[1] == x[3]));
