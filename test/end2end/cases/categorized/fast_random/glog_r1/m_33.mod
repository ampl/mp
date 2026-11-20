var x{1..5} integer >= 0 <= 4;

minimize o: 6*x[1] + x[2] + x[3] + x[4] + 0*x[5];

s.t. c1: (if x[1] == 1 then x[5] else x[4]) > 2;
s.t. c2: (x[4] > 2) ==> (x[3] < 4);
s.t. c3: atleast 1 {i in 1..5} (x[i] >= 3);
