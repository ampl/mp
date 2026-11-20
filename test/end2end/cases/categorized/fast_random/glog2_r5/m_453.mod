var x{1..3} integer >= 0 <= 4;

maximize o: -2*x[1] + 0*x[2] + 3*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) != 2;
s.t. c2: (x[3] == 1) ==> (x[3] == 4) else (x[3] > 1);
s.t. c3: x[3] == 3;
s.t. c4: x[2] == 2;
