var x{1..3} integer >= 0 <= 2;

maximize o: -2000*x[1] - 3000*x[2] - 999*x[3];

s.t. c1: (x[1] == 2) ==> (x[3] == 0) else (x[3] == 2);
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (if x[3] == 0 then x[3] else x[1]) >= 1;
