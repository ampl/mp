\ Model amplgurobi
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x7340a63511dd1bbf
Minimize
  - C0 + C1 - C2 + 0 C3 + 0 C4 + 0 C5 + 0 C6 + 0 C7 + 0 C8 + 0 C9 + 0 C10
   + 0 C11 + 0 C12 + 0 C13 + 0 C14 + 0 C15
Subject To
 R0: C0 + C1 + C2 >= 15
 GC1: C4 = 1 -> 3 C0 - C1 = -19
 GC2: C6 = 1 -> 2 C0 - 2 C2 = -27.99999
 GC3: C8 = 1 -> C1 + 5 C2 = 65
 GC4: C14 = 1 -> C0 = 3.8
 GC5: C3 = 1 -> C1 = 0
 GC6: C12 = 1 -> C0 + C2 <= 2.4
Bounds
 -3 <= C0 <= 7
 C1 <= 10
 C2 <= 11
 C3 <= 1
 C4 <= 1
 C5 = 0
 C6 <= 1
 C7 = 0
 C8 <= 1
 C9 = 0
 C10 <= 1
 C11 = 0
 C12 <= 1
 C13 = 0
 C14 <= 1
 C15 = 1
Generals
 C1 C2 C3 C4 C6 C8 C10 C12 C14
General Constraints
 GC0: C15 = OR ( C3 , C4 , C6 , C10 , C12 , C14 )
\ C10 = 1 + (-1 * (1 + (-1 * C8)))
 GC7: C10 = NL : ( PLUS , -1 , -1 ) ( CONSTANT , 1 , 0 ) 
  ( MULTIPLY , -1 , 0 ) ( CONSTANT , -1 , 2 ) ( PLUS , -1 , 2 ) 
  ( CONSTANT , 1 , 4 ) ( MULTIPLY , -1 , 4 ) ( CONSTANT , -1 , 6 ) 
  ( VARIABLE , C8 , 6 )
End
