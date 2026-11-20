var x{1..3} integer >= 0 <= 2;

minimize o: -4*x[1] - 5*x[2] + 4*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) == 3;
s.t. c2: (if x[1] >= 2 then x[3] else x[1]) >= 0;
s.t. c3: not ((x[1] < 1) and (x[2] == 0));
s.t. c4: alldiff(x[3], x[2], x[1]);
