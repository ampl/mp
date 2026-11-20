var x{1..6} integer >= 0 <= 3;

maximize o: -3000*x[1] - 1000*x[2] + 51*x[3] - 999*x[4] - 2000*x[5] - 3000*x[6];

s.t. c1: not ((x[2] < 2) and (x[6] >= 2));
s.t. c2: max(x[2], x[5]) - min(x[4], x[5]) != 0;
s.t. c3: atleast 3 {i in 1..6} (x[i] >= 3);
s.t. c4: ((if x[4] == 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
