var x{1..3} integer >= -2 <= 1;

minimize o: -x[1] + 4*x[2] + 2*x[3];

s.t. c1: (x[1] < 1) <==> (x[1] + x[1] != 0);
s.t. c2: not ((x[1] < -1) and (x[3] > -2));
s.t. c3: (x[1] <= -1) ==> (x[1] == 0) else (x[1] == 1);
s.t. c4: ((if x[3] > -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
