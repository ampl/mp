var x{1..5} integer >= -1 <= 2;

maximize o:
   -3*x[1] - 2*x[2] - 5*x[3] + 6*x[4] + 2*x[5]
   - 2*(if x[1] < x[1] then 1 else 0);

s.t. c1: (min(x[1], x[1]) < 1) or (abs(x[1] - x[5]) >= 2);
s.t. c2: (abs(x[1] - x[5]) >= 2) ==> (x[1] < x[1]) else (min(x[1], x[1]) < 1);
s.t. c3: (x[1] < x[1]) or (min(x[1], x[1]) < 1);
s.t. c4:
   (if min(x[1], x[1]) < 1 then x[1] else 0)
   + (if abs(x[1] - x[5]) >= 2 then x[2] else 0)
   > -1;
