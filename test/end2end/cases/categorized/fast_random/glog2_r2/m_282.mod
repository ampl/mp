var x{1..6} integer >= -2 <= 0;
maximize o: (-1)*x[1] + (-5)*x[2] + (-1)*x[3] + (-5)*x[4] + (1)*x[5] + (1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: not ((x[3] > -1) and (x[4] < -1));
s.t. c3: atleast 1 {i in 1..6} (x[i] >= -1);
