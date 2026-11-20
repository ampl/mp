var x{1..3} integer >= -2 <= 1;

maximize o: -3*x[1] + 3*x[2] - 3*x[3];

s.t. c1: ((x[1] > 1) ==> (x[3] < 0)) and ((x[3] < 0) ==> (x[3] >= 1));
s.t. c2: (x[1] <= 0) ==> (x[2] <= 1) else (x[3] == -2);
s.t. c3: numberof 0 in (x[1], x[2], x[3]) >= 0;
