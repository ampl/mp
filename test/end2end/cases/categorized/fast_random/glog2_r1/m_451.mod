var x{1..5} integer >= 0 <= 3;

minimize o: 3*x[1] - 4*x[2] + 2*x[3] - 5*x[4] + 4*x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= 0);
s.t. c2: (if x[4] > 0 then x[5] else x[3]) != 4;
s.t. c3: (x[2] < 3) ==> (x[5] >= 0);
s.t. c4: atleast 1 {i in 1..5} (x[i] >= -1);
