var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] + 2*x[2] + x[3];

s.t. c1: (x[1] < 0) ==> (x[1] > 2) else (x[3] >= 3);
s.t. c2: ((x[3] != 2) ==> (x[1] <= -1)) and ((x[1] <= -1) ==> (x[3] >= 2));
s.t. c3: numberof 2 in (x[1], x[2], x[3]) != 3;
