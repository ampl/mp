\ Model amplgurobidirectmodel
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3.19 C0 + 2.59 C1 + 2.29 C2 + 2.89 C3 + 1.89 C4 + 1.99 C5 + 1.99 C6
   + 2.49 C7
Subject To
 R0: 60 C0 + 8 C1 + 8 C2 + 40 C3 + 15 C4 + 70 C5 + 25 C6 + 60 C7 + RgR0
   = 10000
 R1: 10 C0 + 20 C1 + 15 C2 + 35 C3 + 15 C4 + 15 C5 + 25 C6 + 15 C7 + RgR1
   = 10000
 R2: 15 C0 + 20 C1 + 10 C2 + 10 C3 + 15 C4 + 15 C5 + 15 C6 + 10 C7 + RgR2
   = 10000
 R3: 20 C0 + 10 C2 + 40 C3 + 35 C4 + 30 C5 + 50 C6 + 20 C7 + RgR3 = 10000
Bounds
 C0 <= 100
 C1 <= 100
 C2 <= 100
 C3 <= 100
 C4 <= 100
 C5 <= 100
 C6 <= 100
 C7 <= 100
 RgR0 <= 9300
 RgR1 <= 9300
 RgR2 <= 9300
 RgR3 <= 9300
End
