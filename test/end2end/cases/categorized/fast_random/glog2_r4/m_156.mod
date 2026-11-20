var x{1..3} integer >= 0 <= 4;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: atleast 1 {i in 1..3} (x[i] >= 1);
