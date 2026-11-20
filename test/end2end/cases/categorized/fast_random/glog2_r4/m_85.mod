var x{1..3} integer >= 0 <= 2;

minimize o: -4*x[1] - 4*x[2] + 3*x[3];

s.t. c1: (x[3] > 1) ==> (x[1] > 0) else (x[2] >= 1);
s.t. c2: ((x[1] < 0) ==> (x[3] >= 0)) and ((x[3] >= 0) ==> (x[3] >= 0));
s.t. c3: (x[1] != 0) ==> (x[1] < 1);
s.t. c4: atmost 2 {i in 1..3} (x[i] >= 2);
s.t. c5: max(x[2], x[2]) - min(x[2], x[1]) >= 0;
s.t. c6: ((if x[3] < 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c7: (x[1] <= 0) <==> (x[1] + x[2] <= 2);
