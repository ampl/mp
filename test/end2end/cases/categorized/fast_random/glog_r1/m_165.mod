var x{1..3} integer >= -2 <= 1;

maximize o: 5*x[1] + 5*x[2] - 3*x[3];

s.t. c1: abs(x[3] - x[1]) != 3;
s.t. c2: ((x[1] > 0) ==> (x[3] >= -1)) and ((x[3] >= -1) ==> (x[3] >= 0));
s.t. c3: atmost 1 {i in 1..3} (x[i] >= -2);
s.t. c4: (if x[1] > 1 then x[3] else x[1]) > -1;
