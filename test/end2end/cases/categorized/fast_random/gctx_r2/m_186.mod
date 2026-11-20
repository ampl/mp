var x{1..3} integer >= 0 <= 2;

minimize o: x[1] - x[2] + 4*x[3];

s.t. c1: (max(x[2], x[2]) != 2) ==> (max(x[2], x[2]) != 2);
s.t. c2: (max(x[2], x[2]) != 2) ==> (max(x[2], x[2]) != 2);
s.t. c3: (x[3] < 0) <==> (max(x[2], x[2]) != 2);
