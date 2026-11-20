var x{1..3} integer >= -2 <= 3;

maximize o: -3*x[1] - 5*x[2] - 5*x[3];

s.t. c1: ((x[2] == -2) ==> (x[1] <= 1)) and ((x[1] <= 1) ==> (x[3] >= -2));
s.t. c2: (x[1] >= 0) or (x[2] > 0) or (x[2] >= 2);
s.t. c3: (if x[2] == 0 then x[3] else x[1]) < 2;
s.t. c4: max(x[2], x[2]) - min(x[1], x[3]) <= 3;
