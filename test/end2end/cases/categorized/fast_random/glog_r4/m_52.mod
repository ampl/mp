var x{1..3} integer >= -2 <= 2;

maximize o: -x[1] - x[2] - 2*x[3];

s.t. c1: atmost 3 {i in 1..3} (x[i] >= -2);
