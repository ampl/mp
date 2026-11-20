var x{1..4} integer >= 0 <= 4;

maximize o: -3*x[1] - 5*x[2] + 5*x[3] - 2*x[4];

s.t. c1: ((x[2] == 4) ==> (x[2] != 1)) and ((x[2] != 1) ==> (x[2] >= 4));
s.t. c2: abs(x[4] - x[4]) <= 1;
s.t. c3: alldiff(x[2], x[1]);
