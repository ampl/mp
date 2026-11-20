var x{1..5} integer >= -2 <= 3;

maximize o: -4*x[1] - 5*x[2] - x[3] - 2*x[4] - 3*x[5];

s.t. c1: (if x[4] == 3 then x[4] else x[3]) == 1;
s.t. c2: ((if x[4] <= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[4] < -2) ==> (x[3] < 0);
