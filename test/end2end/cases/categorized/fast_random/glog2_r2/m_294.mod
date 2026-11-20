var x{1..3} integer >= -3 <= 1;

minimize o: 5*x[1] + 2*x[2] - x[3];

s.t. c1: numberof -2 in (x[1], x[2], x[3]) <= 2;
s.t. c2: (x[2] >= -2) and (x[1] < -1) ==> (x[2] >= -2);
s.t. c3: alldiff(x[3], x[2]);
