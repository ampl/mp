var x{1..5} integer >= -2 <= 3;

maximize o: -x[1] - 2*x[2] + x[3] + 0*x[4] - x[5];

s.t. c1: (x[1] == 0) ==> (x[1] <= 1) else (x[4] > -1);
s.t. c2: not ((x[4] == -2) and (x[5] <= -2));
s.t. c3: abs(x[3] - x[5]) >= 1;
