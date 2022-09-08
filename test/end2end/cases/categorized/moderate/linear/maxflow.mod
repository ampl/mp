set nodes;
param orig symbolic in nodes;
param dest symbolic in nodes, <> orig;

set arcs within (nodes diff {dest}) cross (nodes diff {orig});

param cap {arcs} >= 0;
var Flow {(i,j) in arcs} >= 0, <= cap[i,j];

maximize Total_Flow:  sum {(orig,j) in arcs} Flow[orig,j];

subject to Balance {k in nodes diff {orig,dest}}:
   sum {(i,k) in arcs} Flow[i,k] = sum {(k,j) in arcs} Flow[k,j];
