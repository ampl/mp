var x{1..4} integer >= 0 <= 2;

minimize o: 2*x[1] + 0*x[2] + x[3] + 3*x[4];

s.t. c1: (x[1] <= 2) ==> (x[2] <= 1);
s.t. c2: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4]) == 0;
