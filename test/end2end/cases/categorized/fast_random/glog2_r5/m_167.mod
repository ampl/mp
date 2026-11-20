var x{1..4} integer >= 0 <= 3;

maximize o: -999*x[1] - 333*x[2] - 3000*x[3] + 2000*x[4];

s.t. c1: not ((x[3] != 2) and (x[2] != 2));
s.t. c2: (x[4] != 0) ==> (x[3] <= 3) else (x[2] < 1);
s.t. c3: abs(x[1] - x[3]) > 3;
