var x{1..4} integer >= 0 <= 2;

maximize o:
   3*x[1] - 2*x[2] + 5*x[3] + 4*x[4] - 6*(if x[1] + x[3] > 2 then 1 else 0);

s.t. c1: (if x[3] == 2 then x[3] else x[4]) >= 2;
s.t. c2: (x[1] + x[3] > 2) <==> (x[3] == 2);
s.t. c3: not (x[1] >= 1);
