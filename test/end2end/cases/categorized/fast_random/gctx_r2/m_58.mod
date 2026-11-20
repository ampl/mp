var x{1..3} integer >= 0 <= 2;

minimize o: 2*x[1] - x[2] - 4*x[3];

s.t. c1: (if x[1] + x[3] <= 0 then x[1] else x[2]) <= 2;
s.t. c2: (x[3] <= 2) ==> ((x[2] < x[3]) or (x[3] <= 2));
s.t. c3:
   (if x[1] + x[3] <= 0 then x[2] else 0) + (if x[3] <= 2 then x[2] else 0)
   <= 0;
