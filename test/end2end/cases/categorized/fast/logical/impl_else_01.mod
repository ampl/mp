## Issue #193

set NS;
set F2 {NS};
set NC;

var X_F2 {s in NS, F2[s], NS, NC};
var U_F2 {s in NS} = sum {v in F2[s], i in NC} X_F2[s,v,s,i];
var VAL_U_F2 {s in NS} binary;

subject to Contrainte_64 {s in NS}:
   VAL_U_F2[s] = 1  ==>  U_F2[s] >= 0.0001  else  U_F2[s] = 0;


data;

set NS := A B C ;
set F2[A] := A1 A3 A4 ;
set F2[B] := B1 B4 ;
set F2[C] := C3 C4 C5 C7 ;
set NC := 11 22 ;

