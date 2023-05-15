Summary of recent updates to MOSEK for AMPL
===========================================


## 20230515
- *Exponential cones*. MP driver recognizes exponential
  cones from their algebraic representation and passes
  them to Mosek.


## 20230505
- *Updated Mosek library* to version 10.0.43. It includes a 
  number of bug fixes, including a numerical problem that 
  could occur with disjunctive constraints


## 20230424
- *Changes in the MP library*: added variable names support
  and removed spurious starting solution
  

## 20230330
- Removed the message on missing dual solution for MIP.


## 20230321
- Recognition of second-order conic constraints
  from algebraic representations.

- First release in the AMPL bundle.


## 20220420
- First release of mock driver.
