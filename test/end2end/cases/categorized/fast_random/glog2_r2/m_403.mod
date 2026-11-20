var x{1..5} integer >= 0 <= 4;

maximize o: 4*x[1] - 5*x[2] + 5*x[3] + 2*x[4] + 5*x[5];

s.t. c1: ((if x[4] < 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c2: not ((x[5] == 3) and (x[4] == 4));
s.t. c3: (x[3] >= 1) and (x[1] != 3) ==> (x[1] >= 1);
s.t. c4: (x[1] != 4) <==> (x[3] + x[2] <= 2);
s.t. c5: x[1] == 3;
