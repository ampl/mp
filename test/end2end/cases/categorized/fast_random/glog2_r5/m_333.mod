var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (4)*x[2] + (3)*x[3];
s.t. c1: ((x[2] == 1) ==> (x[3] > 1)) and ((x[3] > 1) ==> (x[2] >= 1));
s.t. c2: atleast 2 {i in 1..3} (x[i] >= 0);
s.t. c3: x[1] == 1;
s.t. c4: x[2] == 2;
