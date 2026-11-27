var x{1..3} integer >= 0 <= 4;

minimize o: 5*x[1] - 5*x[2] + 2*x[3];

s.t. c1: ((if x[2] >= 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[2] >= 4) and (x[1] < 2) ==> (x[2] >= 4);
s.t. c3: (x[3] <= 4) ==> (x[3] > 1) else (x[1] < 4);
s.t. c4: max(x[2], x[1]) - min(x[1], x[3]) > 3;
s.t. c5: alldiff(x[2], x[1]);
s.t. c6: ((if x[3] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c7: (x[1] <= 2) <==> (x[3] + x[3] == 3);
