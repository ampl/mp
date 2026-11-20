var x{1..3} integer >= -2 <= 3;
minimize o: (-4)*x[1] + (1)*x[2] + (-5)*x[3];
s.t. c1: atleast 2 {i in 1..3} (x[i] >= -2);
s.t. c2: (x[2] < 0) <==> (x[3] + x[2] >= 1);
s.t. c3: alldiff(x[1], x[2]);
