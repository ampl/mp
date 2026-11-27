var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] + 0*x[2] - x[3] + x[4];

s.t. c1: not ((x[2] >= 0) and (x[4] == 0));
s.t. c2: (if x[3] == 2 then x[2] else x[4]) <= 1;
s.t. c3: ((x[2] != 2) ==> (x[3] != 0)) and ((x[3] != 0) ==> (x[1] >= 2));
