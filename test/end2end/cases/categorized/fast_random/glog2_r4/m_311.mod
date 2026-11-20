var x{1..5} integer >= 0 <= 2;

minimize o: -2*x[1] + 0*x[2] - 2*x[3] + 5*x[4] + 5*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) > 1;
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 2);
s.t. c3: (if x[1] != 3 then x[3] else x[3]) < 3;
