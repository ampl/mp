var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] + x[2] - x[3] - 2*(if x[1] < x[2] then 1 else 0);

s.t. c1: not ((x[2] < 2) and (x[2] < 2));
s.t. c2: not ((x[3] > 3) and (x[2] < 2));
s.t. c3: (x[1] < x[2]) <==> (x[2] < 2);
