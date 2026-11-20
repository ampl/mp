var x{1..5} integer >= 0 <= 2;

minimize o: -5*x[1] - 4*x[2] + x[3] - 5*x[4] - 4*x[5];

s.t. c1: not ((x[1] < 0) and (x[3] != 2));
s.t. c2: (if x[3] != 2 then x[1] else x[4]) <= 0;
s.t. c3: (x[2] >= x[2]) ==> ((x[1] < 0) or (x[3] != 2));
