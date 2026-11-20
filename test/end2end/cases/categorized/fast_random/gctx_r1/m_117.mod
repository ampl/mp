var x{1..4} integer >= 0 <= 2;

maximize o: -3*x[1] + 6*x[2] - 3*x[3] - 4*x[4];

s.t. c1: (x[4] == 0) or (x[2] > x[2]);
s.t. c2: (x[2] > x[2]) ==> (x[4] == 0) else (x[2] > x[2]);
s.t. c3: (if x[1] == x[1] then x[1] else x[3]) > 0;
s.t. c4: (x[2] > x[2]) and ((x[4] == 0) or not (x[2] > x[2]));
