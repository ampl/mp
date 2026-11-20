var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] - 5*x[2] + 4*x[3];

s.t. c1: (x[1] == 2) ==> (x[2] < 1) else (x[1] <= 2);
s.t. c2: (if x[2] != 4 then x[1] else x[1]) != 1;
s.t. c3: (if x[2] <= 4 then x[1] else x[3]) > 2;
s.t. c4: abs(x[3] - x[1]) == 1;
s.t. c5: (x[1] == 2) ==> (x[2] < 1) else (x[1] <= 2);
s.t. c6: abs(x[3] - x[1]) == 1;
