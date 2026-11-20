var x{1..5} integer >= 0 <= 2;

minimize o: -x[1] + 5*x[2] - 4*x[3] + 6*x[4] + 2*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c2: ((x[3] < 0) ==> (x[4] != 1)) and ((x[4] != 1) ==> (x[1] >= 0));
