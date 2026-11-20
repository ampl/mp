var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] - 4*x[2] + 5*x[3];

s.t. c1: (if x[2] < 3 then x[2] else x[3]) > 3;
s.t. c2: x[3] == 1;
s.t. c3: x[2] == 3;
