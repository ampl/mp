var x{1..3} integer >= -2 <= 0;

maximize o: -x[1] - 4*x[2] - 5*x[3];

s.t. c1: not ((x[1] <= 0) and (x[1] == 0));
s.t. c2: ((if x[1] >= -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
s.t. c3: (if x[2] <= -1 then x[3] else x[2]) >= -2;
