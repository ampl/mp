var x{1..5} integer >= -2 <= 1;
minimize o: (-3)*x[1] + (2)*x[2] + (-1)*x[3] + (2)*x[4] + (-5)*x[5];
s.t. c1: (x[1] > -2) and (x[1] >= 0) ==> (x[3] >= -2);
s.t. c2: (x[2] >= 0) or (x[3] <= -1) or (x[4] < 0);
s.t. c3: alldiff(x[2], x[4]);
s.t. c4: (if x[3] > -1 then x[1] else x[3]) != 0;
s.t. c5: x[5] == 0;
