var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] - 3*x[2] - 4*x[3];

s.t. c1: (x[2] > 5) ==> (x[2] < 0);
s.t. c2: max(x[3], x[2]) - min(x[3], x[2]) >= 5;
s.t. c3: abs(x[1] - x[2]) >= 3;
