var x{1..3} integer >= -2 <= 2;

minimize o: 0*x[1] - 5*x[2] - x[3];

s.t. c1: not ((x[3] == -2) and (x[1] <= -2));
s.t. c2: ((if x[2] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c3: (if x[3] < -2 then x[1] else x[2]) <= 0;
