var x{1..6} integer >= -2 <= 1;
maximize o: (-2)*x[1] + (-1)*x[2] + (3)*x[3] + (-1)*x[4] + (4)*x[5] + (-3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: exists {i in 1..6} x[i] == 0;
s.t. c3: atleast 3 {i in 1..6} (x[i] >= 1);
s.t. c4: ((if x[3] == -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c5: ((x[3] != -2) ==> (x[4] > 1)) and ((x[4] > 1) ==> (x[2] >= -2));
s.t. c6: (x[1] != -1) or (x[1] != 1) or (x[1] < -2);
s.t. c7: abs(x[5] - x[5]) < 1;
