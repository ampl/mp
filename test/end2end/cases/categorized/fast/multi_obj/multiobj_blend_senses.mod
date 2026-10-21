# Objectives of different senses blended in one priority level.
# The pass objective is the weighted sum in the 1st objective's sense:
#   maximize -3*o2 + 10*o3 = -6x + 7y - 15,
# so that x=0, y=5, and the 2nd pass objective value is 20.
var x integer >= 0, <= 10;
var y integer >= 0, <= 10;

s.t. c: 3*x + 2*y <= 11;

maximize o1: x + y + 10;
minimize o2: 2*x + y + 5;
maximize o3: y;

suffix objpriority IN;
suffix objweight IN;

let o1.objpriority := 2;
let o2.objpriority := 1;
let o3.objpriority := 1;

let o1.objweight := 1;
let o2.objweight := 3;
let o3.objweight := 10;
