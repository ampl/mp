var x{1..3} integer >= -1 <= 2;

minimize o: 3*x[1] + 0*x[2] - 4*x[3];

s.t. c1: (x[2] >= 1) ==> (x[2] < x[2]) else (x[2] < x[2]);
s.t. c2: (x[2] < x[2]) ==> (x[2] >= 1);
s.t. c3: not ((x[2] >= x[2]) and (x[2] >= x[2]));
