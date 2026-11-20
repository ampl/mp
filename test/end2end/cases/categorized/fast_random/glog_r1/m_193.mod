var x{1..5} integer >= 0 <= 3;

minimize o: -3*x[1] - 4*x[2] + 6*x[3] + 5*x[4] + 3*x[5];

s.t. c1: exists {i in 1..5} x[i] == 3;
s.t. c2: ((x[4] < 0) ==> (x[3] > 3)) and ((x[3] > 3) ==> (x[1] >= 0));
s.t. c3: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c4: (if x[5] >= 3 then x[1] else x[3]) < 2;
