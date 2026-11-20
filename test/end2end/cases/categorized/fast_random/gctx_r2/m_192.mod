var x{1..5} integer >= -1 <= 2;

maximize o: -x[1] + 4*x[2] + 6*x[3] - 4*x[4] - 2*x[5];

s.t. c1: not ((x[5] + x[3] == 1) and (x[4] != -1));
s.t. c2: (if x[4] != -1 then x[1] else x[5]) <= -1;
s.t. c3: (x[5] + x[3] == 1) ==> ((x[5] + x[3] == 1) or (x[4] != -1));
s.t. c4: (if x[4] != -1 then x[5] else x[2]) < 0;
