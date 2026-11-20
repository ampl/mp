var x{1..3} integer >= 0 <= 4;

minimize o: 3*x[1] + 3*x[2] + 2*x[3];

s.t. c1: (x[2] != 2) ==> (x[3] != 2);
s.t. c2: max(x[1], x[3]) - min(x[1], x[1]) < 1;
s.t. c3: x[2] == 4;
