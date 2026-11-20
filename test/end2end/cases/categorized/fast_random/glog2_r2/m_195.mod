var x{1..6} integer >= 0 <= 3;

minimize o: 2*x[1] - 4*x[2] - x[3] + 3*x[4] + 5*x[5] + 4*x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 2;
s.t. c2: alldiff(x[2], x[6], x[1], x[5]);
s.t. c3: (x[5] != 1) ==> (x[5] < 1) else (x[2] == 2);
s.t. c4: ((if x[4] == 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
s.t. c5: ((x[2] <= 2) ==> (x[4] <= 3)) and ((x[4] <= 3) ==> (x[4] >= 2));
