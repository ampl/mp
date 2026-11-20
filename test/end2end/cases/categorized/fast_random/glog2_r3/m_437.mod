var x{1..4} integer >= 0 <= 2;
minimize o: (-1000)*x[1] + (-999)*x[2] + (-1000)*x[3] + (-2000)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[3] <= 0) <==> (x[2] + x[1] != 0);
s.t. c3: (if x[4] == 0 then x[4] else x[2]) >= 0;
