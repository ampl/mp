var x{1..5} integer >= -2 <= 0;

maximize o: 0*x[1] - x[2] + 5*x[3] + 2*x[4] + 0*x[5];

s.t. c1: (if x[2] > -2 then x[4] else x[2]) != -2;
s.t. c2: not ((x[2] <= -1) and (x[2] <= -2));
s.t. c3: (if x[4] >= 0 then x[3] else x[4]) != -2;
s.t. c4: (if x[5] >= -1 then x[5] else x[3]) >= -2;
s.t. c5: x[3] == -1;
