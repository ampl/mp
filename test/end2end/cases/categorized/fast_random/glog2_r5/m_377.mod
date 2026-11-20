var x{1..3} integer >= 0 <= 2;

minimize o: 2000*x[1] - 999*x[2] + 17*x[3];

s.t. c1: (if x[3] >= 0 then x[1] else x[1]) > 2;
s.t. c2: (if x[3] < 0 then x[3] else x[3]) < 0;
s.t. c3: (x[3] > 2) or (x[1] >= 2) or (x[1] > 1);
s.t. c4: atleast 1 {i in 1..3} (x[i] >= 0);
