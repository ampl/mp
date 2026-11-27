var x{1..6} integer >= 0 <= 2;

maximize o: 6*x[1] + 2*x[2] + 6*x[3] - 5*x[4] + x[5] - 2*x[6];

s.t. c1: ((if x[5] != 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: atleast 2 {i in 1..6} (x[i] >= 0);
s.t. c3: (x[4] <= 2) ==> (x[4] != 1) else (x[5] == 2);
s.t. c4: ((x[6] < 2) ==> (x[5] <= 1)) and ((x[5] <= 1) ==> (x[4] >= 2));
s.t. c5: not ((x[3] >= 0) and (x[4] < 0));
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
