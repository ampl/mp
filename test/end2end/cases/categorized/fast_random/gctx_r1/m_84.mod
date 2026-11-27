var x{1..5} integer >= 0 <= 3;

maximize o: 4*x[1] - 5*x[2] - x[3] + 2*x[4] + 0*x[5];

s.t. c1: (if abs(x[2] - x[5]) != 1 then x[4] else x[5]) > 1;
s.t. c2:
   (if abs(x[2] - x[5]) != 1 then x[2] else 0) + (if x[4] == 0 then x[5] else 0)
   >= 1;
