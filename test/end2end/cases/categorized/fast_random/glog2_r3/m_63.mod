var x{1..3} integer >= 0 <= 4;
minimize o: (4)*x[1] + (0)*x[2] + (4)*x[3];
s.t. c1: ((x[3] == 4) ==> (x[3] > 3)) and ((x[3] > 3) ==> (x[3] >= 4));
s.t. c2: alldiff(x[3], x[2]);
s.t. c3: atleast 1 {i in 1..3} (x[i] >= 0);
s.t. c4: (x[1] != 2) <==> (x[1] + x[3] == 4);
s.t. c5: x[3] == 2;
