var x{1..6} integer >= 0 <= 4;

minimize o: -5*x[1] - 2*x[2] - 5*x[3] - x[4] + 4*x[5] + x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: ((x[3] < 1) ==> (x[4] > 1)) and ((x[4] > 1) ==> (x[2] >= 1));
s.t. c3: ((x[4] <= 0) ==> (x[4] <= 4)) and ((x[4] <= 4) ==> (x[3] >= 0));
s.t. c4: not ((x[6] != 0) and (x[2] != 2));
s.t. c5: (x[4] <= 3) <==> (x[6] + x[2] == 1);
s.t. c6: ((if x[2] == 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
