var x{1..5} integer >= -1 <= 2;

maximize o:
   -4*x[1] + x[2] + 2*x[3] + 3*x[4] + 4*x[5]
   - 2*(if max(x[1], x[1]) == 1 then 1 else 0);

s.t. c1:
   (max(x[1], x[1]) == 1)
   and ((max(x[1], x[1]) == 1) or not (max(x[1], x[1]) == 1));
s.t. c2: (if x[2] + x[4] == 0 then x[4] else x[5]) < -1;
