var x{1..5} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c2: ((x[1] == 0) ==> (x[1] != 3)) and ((x[1] != 3) ==> (x[3] >= 0));
