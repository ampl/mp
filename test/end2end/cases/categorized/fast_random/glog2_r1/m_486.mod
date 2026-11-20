var x{1..3} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: atleast 3 {i in 1..3} (x[i] >= 1);
s.t. c2: (if x[3] > 2 then x[1] else x[1]) != 4;
