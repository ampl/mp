var x{1..5} integer >= -2 <= 1;

maximize o: -4*x[1] + 5*x[2] - 5*x[3] + 5*x[4] + 3*x[5];

s.t. c1: (x[4] < -2) ==> (x[1] <= 0);
s.t. c2: not ((x[2] <= 1) and (x[5] < 1));
s.t. c3: ((if x[4] != -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
s.t. c4: (x[5] <= 1) and (x[5] == -1) ==> (x[3] >= 1);
s.t. c5: (x[4] <= 0) and (x[5] <= 0) ==> (x[1] >= 0);
s.t. c6: atleast 1 {i in 1..5} (x[i] >= 0);
s.t. c7: (x[2] == -1) ==> (x[2] == -1);
