var x{1..5} integer >= -1 <= 2;

minimize o:
   0*x[1] + 0*x[2] + 3*x[3] + 5*x[4] + 6*x[5]
   + 6*(if x[3] + x[2] != 1 then 1 else 0);

s.t. c1: (x[3] == 1) ==> ((x[1] > x[1]) or (x[3] == 1));
s.t. c2: not (x[1] > x[1]);
s.t. c3: (min(x[4], x[1]) < 0) <==> (x[3] + x[2] != 1);
