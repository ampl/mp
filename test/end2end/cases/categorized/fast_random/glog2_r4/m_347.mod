var x{1..4} integer >= -2 <= 0;

minimize o: -3000*x[1] + 3000*x[2] - 2000*x[3] + 1000*x[4];

s.t. c1: ((if x[2] == 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[4] <= -1) and (x[1] < -2) ==> (x[2] >= -1);
s.t. c3: (x[4] <= 0) ==> (x[1] <= -1);
