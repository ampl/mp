var x{1..5} integer >= -2 <= 0;

maximize o: x[1] + 2*x[2] + x[3] - 4*x[4] + 0*x[5];

s.t. c1: ((x[3] > 1) ==> (x[3] < -2)) and ((x[3] < -2) ==> (x[3] >= 1));
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) != 0;
