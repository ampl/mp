var x{1..5} integer >= 0 <= 3;

minimize o:
   3*x[1] + 3*x[2] + 2*x[3] - 3*x[4] + x[5] - 5*(if x[3] == x[2] then 1 else 0);

s.t. c1: (min(x[1], x[2]) != 3) ==> ((x[3] == x[2]) or (min(x[1], x[2]) != 3));
s.t. c2: (x[3] == x[2]) ==> (x[3] == x[2]);
