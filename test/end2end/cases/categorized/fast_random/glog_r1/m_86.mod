var x{1..3} integer >= 0 <= 2;

minimize o: x[1] - 4*x[2] + 5*x[3];

s.t. c1: (x[1] > 0) ==> (x[2] > 0) else (x[3] == 1);
s.t. c2: abs(x[3] - x[2]) != 1;
s.t. c3: ((if x[1] < 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
