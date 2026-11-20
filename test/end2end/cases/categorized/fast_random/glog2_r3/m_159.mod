var x{1..5} integer >= 0 <= 4;

minimize o: -3*x[1] + 5*x[2] + 2*x[3] - x[4] + 2*x[5];

s.t. c1: (x[5] > 3) ==> (x[3] == 2) else (x[2] != 2);
s.t. c2: not ((x[3] >= 3) and (x[1] <= 2));
s.t. c3: (x[5] > 3) ==> (x[3] == 2) else (x[2] != 2);
s.t. c4: not ((x[3] >= 3) and (x[1] <= 2));
