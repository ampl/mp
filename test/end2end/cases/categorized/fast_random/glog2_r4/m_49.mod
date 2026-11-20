var x{1..5} integer >= -2 <= 2;

maximize o: -5*x[1] + 4*x[2] + 2*x[3] + 5*x[4] - x[5];

s.t. c1: ((if x[2] <= -2 then 1 else 0) + (if x[5] >= -2 then 1 else 0)) >= 1;
s.t. c2: atmost 3 {i in 1..5} (x[i] >= -2);
s.t. c3: (x[4] <= 0) or (x[3] > -2) or (x[1] > 1);
s.t. c4: ((if x[2] <= -2 then 1 else 0) + (if x[5] >= -2 then 1 else 0)) >= 1;
s.t. c5: (x[4] <= 0) or (x[3] > -2) or (x[1] > 1);
