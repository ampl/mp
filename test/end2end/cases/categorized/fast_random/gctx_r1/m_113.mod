var x{1..5} integer >= -1 <= 2;

maximize o: -3*x[1] + 5*x[2] - 2*x[3] + 6*x[4] - x[5];

s.t. c1: (if x[1] + x[2] > -1 then x[1] else x[4]) <= -1;
s.t. c2: (x[1] != 1) ==> (x[1] != 1);
s.t. c3: (x[1] != 1) or (x[1] <= 2);
