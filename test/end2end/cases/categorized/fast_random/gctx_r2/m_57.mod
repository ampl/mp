var x{1..5} integer >= 0 <= 2;

minimize o: 3*x[1] - 3*x[2] - x[3] + 0*x[4] + 0*x[5];

s.t. c1: (x[4] + x[1] < 1) ==> (x[2] >= 1) else (x[2] >= 1);
s.t. c2: (x[3] == x[4]) and ((x[3] == x[4]) or not (x[2] >= 1));
s.t. c3: (x[2] >= 1) or (x[2] >= 1);
s.t. c4: not ((x[3] == x[4]) and (x[4] + x[1] < 1));
