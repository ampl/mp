var x{1..5} integer >= 0 <= 3;

maximize o: 5*x[1] + 3*x[2] + x[3] - 3*x[4] - 2*x[5];

s.t. c1: ((x[1] < 0) ==> (x[2] <= 1)) and ((x[2] <= 1) ==> (x[2] >= 0));
s.t. c2: abs(x[3] - x[3]) < 3;
s.t. c3: not ((x[2] != 3) and (x[3] < 2));
s.t. c4: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) < 0;
