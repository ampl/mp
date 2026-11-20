var x{1..5} integer >= 0 <= 4;

minimize o: -x[1] + 5*x[2] + 0*x[3] + 3*x[4] + 4*x[5];

s.t. c1: (x[4] == 3) ==> (x[4] == 4) else (x[5] <= 1);
s.t. c2: max(x[1], x[5]) - min(x[2], x[1]) >= 0;
s.t. c3: not ((x[2] > 3) and (x[4] == 2));
s.t. c4: x[4] == 0;
s.t. c5: x[3] == 0;
