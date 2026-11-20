var x{1..3} integer >= 0 <= 2;

maximize o: 5*x[1] - x[2] + 0*x[3];

s.t. c1: ((if x[1] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: atmost 1 {i in 1..3} (x[i] >= 2);
