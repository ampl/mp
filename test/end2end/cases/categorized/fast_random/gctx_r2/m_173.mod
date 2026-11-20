var x{1..5} integer >= 0 <= 3;

maximize o: -x[1] + 5*x[2] - 3*x[3] + x[4] + 5*x[5];

s.t. c1: (if x[2] == 0 then x[4] else x[1]) >= 3;
s.t. c2: (x[1] == 2) and ((x[2] == 0) or not (x[1] == 2));
s.t. c3: (if x[2] == 0 then x[2] else 0) + (if x[2] == 0 then x[5] else 0) > 3;
