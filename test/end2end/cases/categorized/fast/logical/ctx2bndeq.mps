* FICO Xpress v9.8.1, Hyper, written 16:28:18, May 18, 2026
* Written from user model (presolvestate=0x680061)
NAME          
ROWS
 N  __OBJ___        
 E  R1              
 E  R2              
 E  R3              
 E  R4              
 E  R5              
 L  R6              
 G  R7              
 G  R8              
COLUMNS
    z                 __OBJ___          -1
    z                 R1                3
    z                 R2                2
    z                 R4                1
    z                 R6                1
    z                 R7                1
    x                 __OBJ___          1
    x                 R1                -1
    x                 R3                1
    x                 R5                1
    x                 R7                1
    y                 __OBJ___          -1
    y                 R2                -2
    y                 R3                5
    y                 R6                1
    y                 R7                1
    LC1_flatv1_       __OBJ___          0
    LC1_flatv2_       __OBJ___          0
    LC1_flatv3_       __OBJ___          0
    LC1_flatv4_       __OBJ___          0
    LC1_flatv5_       __OBJ___          0
    LC1_flatv6_       __OBJ___          0
    LC1_flatv7_       __OBJ___          0
    LC1_flatv8_       R8                -1
    LC1_flatv9_       __OBJ___          0
    LC1_flatv10_      __OBJ___          0
    LC1_flatv11_      __OBJ___          0
    LC1_flatv12_      __OBJ___          0
    LC1_flatv13_      __OBJ___          0
    =                 __OBJ___          0
    =         R8        = ( LC1_flatv6_ * -1 + 1 ) * -1 + 1
RHS
    RHS00001          R1                -19
    RHS00001          R2                -27.99999
    RHS00001          R3                65
    RHS00001          R4                3.8
    RHS00001          R6                2.4
    RHS00001          R7                15
BOUNDS
 UP BND00001          z                 7
 LO BND00001          z                 -3
 UI BND00001          x                 10
 UI BND00001          y                 11
 BV BND00001          LC1_flatv1_     
 BV BND00001          LC1_flatv2_     
 FX BND00001          LC1_flatv3_       0
 BV BND00001          LC1_flatv4_     
 FX BND00001          LC1_flatv5_       0
 BV BND00001          LC1_flatv6_     
 FX BND00001          LC1_flatv7_       0
 BV BND00001          LC1_flatv8_     
 FX BND00001          LC1_flatv9_       0
 BV BND00001          LC1_flatv10_    
 FX BND00001          LC1_flatv11_      0
 BV BND00001          LC1_flatv12_    
 FX BND00001          LC1_flatv13_      1
 FX BND00001          =                 1
INDICATORS
 IF R1                LC1_flatv2_       1
 IF R2                LC1_flatv4_       1
 IF R3                LC1_flatv6_       1
 IF R4                LC1_flatv12_      1
 IF R5                LC1_flatv1_       1
 IF R6                LC1_flatv10_      1
GENCONS
 OR GC1             
    LC1_flatv13_    
    LC1_flatv1_     
    LC1_flatv2_     
    LC1_flatv4_     
    LC1_flatv8_     
    LC1_flatv10_    
    LC1_flatv12_    
SLPDATA
 DR LC1_flatv8_  R8                    
ENDATA
