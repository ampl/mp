var x{1..3} integer >= 0 <= 2;

minimize o: 0*x[1] + 5*x[2] + 2*x[3];

s.t. c1: numberof 0 in (x[1], x[2], x[3]) >= 1;
s.t. c2: (x[1] > 1) ==> (x[1] >= 1) else (x[2] >= 2);
s.t. c3: ((x[3] == 2) ==> (x[1] >= 2)) and ((x[1] >= 2) ==> (x[1] >= 2));
s.t. c4: exists {i in 1..3} x[i] == 1;
