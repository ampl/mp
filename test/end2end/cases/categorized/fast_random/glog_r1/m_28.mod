var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (6)*x[2] + (2)*x[3];
s.t. c1: (x[2] < 2) ==> (x[3] >= 2) else (x[3] < 2);
s.t. c2: (x[2] <= 0) ==> (x[1] < 2);
s.t. c3: alldiff(x[1], x[3], x[2]);
s.t. c4: atleast 1 {i in 1..3} (x[i] >= 0);
