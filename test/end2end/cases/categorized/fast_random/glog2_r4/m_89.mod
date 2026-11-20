var x{1..5} integer >= -2 <= 0;

maximize o: 2*x[1] - 2*x[2] + 6*x[3] + 6*x[4] - 5*x[5];

s.t. c1: max(x[3], x[3]) - min(x[1], x[5]) != -2;
s.t. c2: ((if x[5] >= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: abs(x[4] - x[1]) <= 0;
s.t. c4: (x[5] >= -2) and (x[1] <= -2) ==> (x[2] >= -2);
s.t. c5: max(x[3], x[3]) - min(x[1], x[5]) != -2;
s.t. c6: (x[5] >= -2) and (x[1] <= -2) ==> (x[2] >= -2);
