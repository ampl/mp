var x{1..6} integer >= 0 <= 2;

maximize o: -4*x[1] - 5*x[2] + 0*x[3] - 4*x[4] - 2*x[5] + 2*x[6];

s.t. c1: (if x[3] < 2 then x[6] else x[3]) <= 2;
s.t. c2: not ((x[4] == 1) and (x[6] <= 1));
s.t. c3: abs(x[2] - x[2]) != 1;
s.t. c4: (if x[3] < 2 then x[6] else x[3]) <= 2;
s.t. c5: abs(x[2] - x[2]) != 1;
