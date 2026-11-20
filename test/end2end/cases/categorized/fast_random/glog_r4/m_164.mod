var x{1..4} integer >= -2 <= 2;

minimize o: -x[1] - 5*x[2] + 6*x[3] - 4*x[4];

s.t. c1: ((if x[1] < 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[3] == -1) ==> (x[2] <= 2);
s.t. c3: (x[1] > 2) <==> (x[1] + x[4] > -1);
