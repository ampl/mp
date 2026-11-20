var x{1..3} integer >= -2 <= 0;

maximize o: 4*x[1] + 4*x[2] - x[3];

s.t. c1: (x[1] != -1) <==> (x[3] + x[2] != -1);
s.t. c2: atleast 2 {i in 1..3} (x[i] >= -2);
s.t. c3: (x[3] < -2) <==> (x[1] + x[3] < -1);
