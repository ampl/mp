var x{1..5} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (-1)*x[2] + (-5)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: not ((x[5] + x[3] > 0) and (x[1] >= x[5]));
s.t. c2: (if x[5] + x[3] > 0 then x[1] else x[5]) >= 1;
s.t. c3: (x[1] >= x[5]) ==> (x[5] + x[3] > 0);
s.t. c4: (if x[5] + x[3] > 0 then x[2] else 0) + (if min(x[4], x[2]) > 2 then x[4] else 0) >= 1;
