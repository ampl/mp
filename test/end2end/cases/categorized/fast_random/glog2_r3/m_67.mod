var x{1..4} integer >= 0 <= 2;
minimize o: (51)*x[1] + (51)*x[2] + (-3000)*x[3] + (3000)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
