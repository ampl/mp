var x{1..3} integer >= 0 <= 4;

maximize o: -2000*x[1] - 333*x[2] + 3000*x[3];

s.t. c1: ((x[3] <= 0) ==> (x[2] <= 2)) and ((x[2] <= 2) ==> (x[3] >= 0));
s.t. c2: not ((x[3] >= 0) and (x[2] != 3));
s.t. c3: numberof 1 in (x[1], x[2], x[3]) <= 1;
s.t. c4: max(x[1], x[3]) - min(x[3], x[1]) == 4;
