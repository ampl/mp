var x{1..3} integer >= 0 <= 4;

maximize o: 2*x[1] + x[2] - 2*x[3];

s.t. c1: max(x[1], x[2]) - min(x[1], x[1]) <= 5;
s.t. c2: numberof 4 in (x[1], x[2], x[3]) <= 0;
s.t. c3: numberof 4 in (x[1], x[2], x[3]) <= 0;
s.t. c4: ((x[1] == 4) ==> (x[2] <= 5)) and ((x[2] <= 5) ==> (x[2] >= 4));
