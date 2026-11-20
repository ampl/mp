var x{1..3} integer >= 0 <= 3;

maximize o: -x[1] - 2*x[2] + 5*x[3];

s.t. c1: ((if x[2] <= 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[3] >= 1) and (x[2] >= 0) ==> (x[3] >= 1);
s.t. c3: x[3] == 2;
s.t. c4: x[2] == 3;
