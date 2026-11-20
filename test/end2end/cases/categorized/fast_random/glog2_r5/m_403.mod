var x{1..5} integer >= -2 <= 0;

maximize o: x[1] + 5*x[2] + 0*x[3] + 6*x[4] - x[5];

s.t. c1: ((if x[5] != 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[2] <= -1) or (x[5] < 0) or (x[3] != -1);
s.t. c3: (x[5] != 0) ==> (x[1] != -1) else (x[5] == -2);
s.t. c4: alldiff(x[3], x[4]);
s.t. c5: x[5] == -2;
