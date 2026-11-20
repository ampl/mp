var x{1..3} integer >= -2 <= 2;

minimize o: -3*x[1] - x[2] + 0*x[3];

s.t. c1: ((x[2] < 2) ==> (x[2] <= 0)) and ((x[2] <= 0) ==> (x[1] >= 2));
s.t. c2: (if x[3] == -2 then x[2] else x[1]) <= -1;
s.t. c3: ((x[2] < 2) ==> (x[2] <= 0)) and ((x[2] <= 0) ==> (x[1] >= 2));
s.t. c4: (if x[3] == -2 then x[2] else x[1]) <= -1;
