var x{1..3} integer >= -2 <= 1;

maximize o: -x[1] + 3*x[2] - x[3];

s.t. c1: exists {i in 1..3} x[i] == -1;
s.t. c2: (if x[3] != 0 then x[2] else x[2]) > 1;
s.t. c3: (x[1] == -2) and (x[3] == 0) ==> (x[1] >= -2);
