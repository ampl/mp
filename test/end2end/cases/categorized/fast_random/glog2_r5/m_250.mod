var x{1..5} integer >= -2 <= 1;

maximize o: 2*x[1] + 2*x[2] + 6*x[3] - 2*x[4] + x[5];

s.t. c1: ((x[1] < -2) ==> (x[1] != -1)) and ((x[1] != -1) ==> (x[1] >= -2));
s.t. c2: atmost 3 {i in 1..5} (x[i] >= -1);
s.t. c3: ((if x[1] != -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c4: (x[2] >= -2) <==> (x[4] + x[3] == -1);
