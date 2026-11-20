var x{1..4} integer >= 0 <= 3;

maximize o: -x[1] + 6*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[4] > x[4]) ==> ((x[4] > x[4]) or (x[4] > x[4]));
s.t. c2: (if x[4] == 1 then x[4] else x[1]) > 3;
