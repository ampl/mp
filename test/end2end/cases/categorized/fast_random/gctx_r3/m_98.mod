var x{1..5} integer >= 0 <= 2;

maximize o:
   0*x[1] - 3*x[2] + 6*x[3] - 3*x[4] + 0*x[5]
   - 5*(if abs(x[4] - x[1]) != 1 then 1 else 0);

s.t. c1: (x[3] <= x[5]) <==> (x[3] <= x[5]);
s.t. c2: (abs(x[4] - x[1]) != 1) <==> (x[3] >= 1);
s.t. c3: not ((x[3] <= x[5]) and (x[3] >= 1));
