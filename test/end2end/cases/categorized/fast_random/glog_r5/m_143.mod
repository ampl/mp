var x{1..3} integer >= 0 <= 4;

maximize o: 3*x[1] + 4*x[2] + x[3];

s.t. c1: (x[1] < 1) or (x[1] != 0) or (x[1] == 0);
s.t. c2: numberof 3 in (x[1], x[2], x[3]) == 2;
s.t. c3: ((x[3] != 3) ==> (x[2] > 2)) and ((x[2] > 2) ==> (x[1] >= 3));
s.t. c4: max(x[3], x[2]) - min(x[3], x[3]) >= 2;
