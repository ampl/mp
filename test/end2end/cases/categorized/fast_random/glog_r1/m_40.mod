var x{1..5} integer >= -2 <= 0;

minimize o: 6*x[1] - 5*x[2] - 3*x[3] + 0*x[4] + x[5];

s.t. c1: (if x[5] >= 0 then x[4] else x[5]) != -2;
s.t. c2: ((x[2] < 0) ==> (x[1] != 0)) and ((x[1] != 0) ==> (x[3] >= 0));
