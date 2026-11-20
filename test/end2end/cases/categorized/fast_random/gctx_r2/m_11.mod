var x{1..4} integer >= 0 <= 2;

maximize o:
   6*x[1] + x[2] - x[3] - 2*x[4] - 2*(if x[1] + x[3] < 1 then 1 else 0);

s.t. c1: (x[2] > 0) <==> (x[2] > 0);
s.t. c2: (x[2] > 0) ==> (x[1] + x[3] < 1) else (x[1] + x[3] < 1);
s.t. c3: (x[2] > 0) ==> ((x[2] > 0) or (x[1] + x[3] < 1));
s.t. c4: (if x[1] + x[3] < 1 then x[3] else x[4]) > 0;
