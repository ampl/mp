var x{1..5} integer >= 0 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c2: not ((x[2] >= 0) and (x[5] == 0));
s.t. c3: (x[3] != 2) ==> (x[2] != 1);
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) != 2;
