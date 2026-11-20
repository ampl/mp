var x{1..4} integer >= 0 <= 2;

maximize o: -x[1] - x[2] - 3*x[3] + 0*x[4];

s.t. c1: (if abs(x[3] - x[3]) > 0 then x[4] else x[2]) < 2;
s.t. c2: (max(x[4], x[1]) <= 2) ==> (max(x[4], x[1]) <= 2);
