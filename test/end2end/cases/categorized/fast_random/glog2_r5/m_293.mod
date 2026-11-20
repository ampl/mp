var x{1..3} integer >= -2 <= 0;

minimize o: 2*x[1] - x[2] + x[3];

s.t. c1: atleast 2 {i in 1..3} (x[i] >= -2);
s.t. c2: atmost 3 {i in 1..3} (x[i] >= 0);
s.t. c3: (x[1] != 0) ==> (x[2] < 0);
s.t. c4: not ((x[1] <= 0) and (x[3] >= -2));
s.t. c5: x[3] == -1;
s.t. c6: x[1] == -1;
