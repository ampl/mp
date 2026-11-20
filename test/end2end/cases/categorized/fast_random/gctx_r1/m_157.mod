var x{1..3} integer >= 0 <= 2;

minimize o: 6*x[1] + 2*x[2] - 2*x[3] + 2*(if min(x[1], x[3]) < 1 then 1 else 0);

s.t. c1: (min(x[1], x[3]) < 1) ==> ((x[2] < 1) or (min(x[3], x[2]) != 2));
s.t. c2: (min(x[3], x[2]) != 2) ==> (x[2] < 1) else (abs(x[1] - x[1]) <= 2);
s.t. c3: not (min(x[3], x[2]) != 2);
