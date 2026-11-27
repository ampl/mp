var x{1..3} integer >= -2 <= 2;

maximize o: 4*x[1] + x[2] - 5*x[3];

s.t. c1: (if x[2] <= 1 then x[3] else x[1]) >= 0;
s.t. c2: (x[1] <= 2) and (x[2] < 0) ==> (x[2] >= 2);
s.t. c3: (x[2] >= -1) or (x[3] != 1) or (x[1] > 2);
