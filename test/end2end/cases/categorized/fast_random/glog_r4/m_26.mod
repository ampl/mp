var x{1..3} integer >= 0 <= 2;

minimize o: -2*x[1] - 5*x[2] + 3*x[3];

s.t. c1: (x[3] >= 0) or (x[3] >= 0) or (x[2] >= 2);
s.t. c2: ((if x[3] == 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c3: (if x[1] < 2 then x[3] else x[3]) > 1;
s.t. c4: (if x[2] != 1 then x[2] else x[3]) >= 0;
