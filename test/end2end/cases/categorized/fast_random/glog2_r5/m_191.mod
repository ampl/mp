var x{1..5} integer >= 0 <= 3;

minimize o: 3*x[1] + 6*x[2] - 3*x[3] - x[4] - 5*x[5];

s.t. c1: ((if x[3] <= -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[2], x[5], x[3], x[1]);
s.t. c3: ((x[2] <= 3) ==> (x[1] != -1)) and ((x[1] != -1) ==> (x[5] >= 3));
s.t. c4: ((if x[5] >= -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
