var x{1..4} integer >= 0 <= 4;
minimize o: (-1000)*x[1] + (17)*x[2] + (34)*x[3] + (-333)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
