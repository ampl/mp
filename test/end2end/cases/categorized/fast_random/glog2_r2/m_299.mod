var x{1..5} integer >= 0 <= 2;

minimize o: x[1] - 3*x[2] + 5*x[3] + 6*x[4] - x[5];

s.t. c1: (x[3] >= 0) <==> (x[3] + x[5] != 1);
s.t. c2: (if x[1] == 2 then x[5] else x[5]) != 1;
s.t. c3: (x[3] >= 0) <==> (x[3] + x[5] != 1);
s.t. c4: (if x[1] == 2 then x[5] else x[5]) != 1;
