\ Model amplgurobi
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x75369a04bc931c59
Minimize
  3 C0 - 2 C1
Subject To
 R0: 3 C0 - 2 C1 = 13
 R1: 2 C2 - 4 C3 <= -34
Bounds
 -2 <= C0 <= 5
 -12 <= C1 <= 3
 C2 free
 C3 free
General Constraints
 GC0: C2 = ABS ( C0 )
 GC1: C3 = ABS ( C1 )
End
