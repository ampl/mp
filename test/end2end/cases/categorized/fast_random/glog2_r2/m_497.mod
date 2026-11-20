var x{1..5} integer >= 0 <= 4;
maximize o: (3000)*x[1] + (3000)*x[2] + (2000)*x[3] + (3000)*x[4] + (-333)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
