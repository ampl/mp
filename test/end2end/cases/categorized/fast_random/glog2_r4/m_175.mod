var x{1..3} integer >= 0 <= 3;

maximize o: 6*x[1] + x[2] + 0*x[3];

s.t. c1: ((x[1] > 2) ==> (x[1] >= 1)) and ((x[1] >= 1) ==> (x[2] >= 2));
s.t. c2: ((if x[2] < 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c3: numberof 2 in (x[1], x[2], x[3]) == 1;
s.t. c4: (x[3] <= 1) or (x[2] != 2) or (x[2] == 1);
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
