var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] - 2*x[2] - 3*x[3];

s.t. c1: numberof 1 in (x[1], x[2], x[3]) != 3;
s.t. c2: not ((x[2] > 0) and (x[2] <= 3));
s.t. c3: atmost 3 {i in 1..3} (x[i] >= 0);
s.t. c4: ((x[3] < 3) ==> (x[3] > 1)) and ((x[3] > 1) ==> (x[1] >= 3));
