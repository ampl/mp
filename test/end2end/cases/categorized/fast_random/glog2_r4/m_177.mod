var x{1..4} integer >= 0 <= 2;
maximize o: (-999)*x[1] + (3000)*x[2] + (2000)*x[3] + (-3000)*x[4];
s.t. c1: ((if x[3] >= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
