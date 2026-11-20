var x{1..3} integer >= -2 <= 1;

minimize o: 3000*x[1] - 999*x[2] - 666*x[3];

s.t. c1: (x[1] < 1) and (x[1] < 0) ==> (x[2] >= 1);
s.t. c2: numberof -2 in (x[1], x[2], x[3]) == 3;
