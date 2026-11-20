var x{1..7} integer >= 0 <= 1;

maximize o: 3*x[1] + 2*x[2] - 3*x[3] - 5*x[4] - x[5] + 3*x[6] + 5*x[7];

s.t. c1: not ((x[2] != 1) and (x[5] != 1));
s.t. c2: (if x[2] > 0 then x[3] else x[5]) != 0;
