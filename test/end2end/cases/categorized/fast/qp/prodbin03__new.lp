\ Model amplgurobi
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x18924b7838510576
Minimize
  b[1] - b[2] + b[3]
Subject To
 C1: C1_auxv2_ = 1
 C1_PHo: - b[2] - C1_auxv1_ = -1
Bounds
 b[1] = 1
 b[2] <= 1
 b[3] = 1
 C1_auxv1_ = 1
 C1_auxv2_ = 1
Generals
 b[2]
End
