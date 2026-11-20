var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] - 2*x[3] - x[4] - 3*(if x[3] > x[2] then 1 else 0);

s.t. c1: (min(x[1], x[3]) == 3) and ((x[4] + x[1] == 1) or not (x[1] != x[4]));
s.t. c2: (if x[3] > x[2] then x[4] else x[2]) < 2;
s.t. c3: (x[4] + x[1] == 1) <==> (x[4] + x[1] == 1);
