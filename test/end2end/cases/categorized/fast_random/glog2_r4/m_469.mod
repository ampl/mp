var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] + 4*x[2] + x[3] + 6*x[4];

s.t. c1: max(x[3], x[3]) - min(x[4], x[3]) < 2;
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c3: (x[3] >= 2) and (x[1] < 0) ==> (x[1] >= 2);
s.t. c4: exists {i in 1..4} x[i] == 3;
s.t. c5: max(x[3], x[3]) - min(x[4], x[3]) < 2;
s.t. c6: exists {i in 1..4} x[i] == 3;
