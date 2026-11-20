var x{1..3} integer >= -2 <= 2;

maximize o: 0*x[1] - 4*x[2] - 5*x[3];

s.t. c1: (x[2] < 1) and (x[2] == 1) ==> (x[3] >= 1);
s.t. c2: (if x[3] >= 0 then x[3] else x[3]) < 1;
s.t. c3: not ((x[1] > 2) and (x[3] != 1));
s.t. c4: exists {i in 1..3} x[i] == 0;
