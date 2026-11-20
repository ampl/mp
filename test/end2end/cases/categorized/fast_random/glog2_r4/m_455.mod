var x{1..4} integer >= -2 <= 1;

minimize o: x[1] - x[2] + 4*x[3] + 0*x[4];

s.t. c1: abs(x[1] - x[3]) < 0;
s.t. c2: (x[2] >= 0) <==> (x[1] + x[1] <= 0);
s.t. c3: (x[3] != -2) <==> (x[2] + x[2] >= -1);
s.t. c4: max(x[2], x[2]) - min(x[3], x[4]) != -2;
s.t. c5: ((if x[4] < -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c6: ((if x[1] > -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c7: (x[4] <= 0) or (x[4] >= -1) or (x[1] >= -2);
