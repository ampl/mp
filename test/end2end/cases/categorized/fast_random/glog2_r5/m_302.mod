var x{1..5} integer >= 0 <= 3;

minimize o: -2*x[1] - 3*x[2] - 4*x[3] + 0*x[4] + 4*x[5];

s.t. c1: ((x[1] < 0) ==> (x[3] <= 1)) and ((x[3] <= 1) ==> (x[1] >= 0));
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) < 1;
s.t. c3: (x[2] >= 3) or (x[5] < 3) or (x[5] == 0);
s.t. c4: (if x[5] <= 2 then x[5] else x[3]) <= 1;
