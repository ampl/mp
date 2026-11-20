var x{1..3} integer >= 0 <= 3;

maximize o: -2*x[1] + 0*x[2] + x[3] + 3*(if x[2] <= 0 then 1 else 0);

s.t. c1: (x[2] <= 0) or (abs(x[3] - x[1]) < 0);
s.t. c2: not ((abs(x[3] - x[1]) < 0) and (abs(x[3] - x[1]) < 0));
s.t. c3:
   (if abs(x[3] - x[1]) < 0 then 1 else 0) + (if x[2] <= 0 then 1 else 0)
   + (if abs(x[3] - x[1]) < 0 then 1 else 0)
   >= 2;
s.t. c4: (x[2] <= 0) <==> (abs(x[3] - x[1]) < 0);
