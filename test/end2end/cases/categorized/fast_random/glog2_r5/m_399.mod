var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] - 3*x[2] + 6*x[3];

s.t. c1: alldiff(x[2], x[3]);
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 3);
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 3);
s.t. c4: alldiff(x[2], x[3]);
s.t. c5: atmost 2 {i in 1..3} (x[i] >= 3);
