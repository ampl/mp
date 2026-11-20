var x{1..3} integer >= 0 <= 3;

maximize o: 4*x[1] - 4*x[2] - 4*x[3];

s.t. c1:
   (abs(x[2] - x[2]) > 2)
   and ((abs(x[3] - x[1]) == 2) or not (abs(x[2] - x[2]) > 2));
s.t. c2: (if abs(x[2] - x[2]) > 2 then x[2] else x[2]) == 1;
s.t. c3:
   (if abs(x[2] - x[2]) > 2 then x[2] else 0)
   + (if abs(x[3] - x[1]) == 2 then x[2] else 0)
   > 1;
