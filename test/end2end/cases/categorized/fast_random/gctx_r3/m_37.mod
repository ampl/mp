var x{1..4} integer >= 0 <= 3;

maximize o: 3*x[1] - 2*x[2] + 5*x[3] + 3*x[4];

s.t. c1:
   (if x[2] >= 0 then x[3] else 0) + (if max(x[2], x[3]) < 3 then x[4] else 0)
   < 2;
s.t. c2: (x[2] >= 0) and ((max(x[2], x[3]) < 3) or not (x[2] >= 0));
