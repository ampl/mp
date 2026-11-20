var x{1..5} integer >= 0 <= 2;

maximize o: -3000*x[1] - 1000*x[2] + 2000*x[3] + 17*x[4] - 2000*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) <= 2;
s.t. c2: abs(x[3] - x[5]) >= 1;
s.t. c3: (x[2] != 0) ==> (x[4] != 0) else (x[5] > 0);
