var x{1..6} integer >= -2 <= 2;

maximize o: 3*x[1] - 5*x[2] - 3*x[3] - 4*x[4] + 3*x[5] + x[6];

s.t. c1: atleast 3 {i in 1..6} (x[i] >= -2);
s.t. c2: ((if x[6] >= -1 then 1 else 0) + (if x[6] >= -1 then 1 else 0)) >= 1;
s.t. c3: abs(x[2] - x[5]) < 0;
