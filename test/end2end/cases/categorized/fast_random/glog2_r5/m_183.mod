var x{1..5} integer >= 0 <= 2;

maximize o: -2*x[1] + 0*x[2] - 5*x[3] + 2*x[4] - 3*x[5];

s.t. c1: ((x[2] == 0) ==> (x[3] > 2)) and ((x[3] > 2) ==> (x[2] >= 0));
s.t. c2: (if x[3] != 2 then x[2] else x[1]) <= 0;
s.t. c3: x[5] == 1;
s.t. c4: x[2] == 1;
