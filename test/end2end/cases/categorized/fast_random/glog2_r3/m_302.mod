var x{1..5} integer >= 0 <= 4;

minimize o: -2*x[1] - 2*x[2] + x[3] - 4*x[4] + 5*x[5];

s.t. c1: (if x[2] == 0 then x[2] else x[2]) > 2;
s.t. c2: (x[3] < 0) ==> (x[3] == 1) else (x[3] != 4);
