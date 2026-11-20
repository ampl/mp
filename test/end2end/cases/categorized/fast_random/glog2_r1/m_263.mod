var x{1..6} integer >= -2 <= 2;

minimize o: 6*x[1] - 2*x[2] + 0*x[3] - 3*x[4] + 0*x[5] - 4*x[6];

s.t. c1: abs(x[3] - x[6]) <= 0;
s.t. c2: (x[5] != -1) or (x[3] <= 2) or (x[4] > -2);
s.t. c3: ((x[3] >= 0) ==> (x[6] < 0)) and ((x[6] < 0) ==> (x[4] >= 0));
s.t. c4: ((if x[4] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c5: x[3] == 2;
s.t. c6: x[2] == -1;
