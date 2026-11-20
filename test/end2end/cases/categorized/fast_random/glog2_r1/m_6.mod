var x{1..3} integer >= 0 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: max(x[3], x[2]) - min(x[2], x[3]) == 2;
s.t. c2: atleast 1 {i in 1..3} (x[i] >= 1);
s.t. c3: ((if x[2] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c4: ((x[1] != 0) ==> (x[3] != 1)) and ((x[3] != 1) ==> (x[1] >= 0));
