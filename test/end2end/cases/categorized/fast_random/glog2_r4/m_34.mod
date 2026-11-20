var x{1..5} integer >= -2 <= 2;

maximize o: 3*x[1] - x[2] + 0*x[3] - 5*x[4] + 2*x[5];

s.t. c1: ((if x[5] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: (if x[5] == 2 then x[4] else x[3]) > 1;
s.t. c3: ((x[5] <= -2) ==> (x[1] < -1)) and ((x[1] < -1) ==> (x[1] >= -2));
