var x{1..5} integer >= -2 <= 1;

minimize o: 2*x[1] - x[2] + 4*x[3] + 6*x[4] - x[5];

s.t. c1: (x[3] > -1) or (x[3] == 1) or (x[4] <= 1);
s.t. c2: (if x[5] < -2 then x[1] else x[4]) >= -2;
s.t. c3: (x[3] > -1) or (x[3] == 1) or (x[4] <= 1);
s.t. c4: (if x[5] < -2 then x[1] else x[4]) >= -2;
