var x{1..6} integer >= -2 <= 1;

minimize o: -5*x[1] + x[2] - 5*x[3] - 3*x[4] + 4*x[5] + 4*x[6];

s.t. c1: (x[1] <= 1) or (x[4] < 1) or (x[4] < 0);
s.t. c2: ((if x[3] > -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c3: (x[1] >= 0) ==> (x[6] == -2);
s.t. c4: (x[2] >= 0) <==> (x[1] + x[6] <= -1);
s.t. c5: x[5] == 0;
