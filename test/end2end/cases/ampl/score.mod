# SCORE: a 2-class credit scoring model

# Find the set of weights that best classifies known
# "good" and "bad" credit risks, in the sense of minimizing a
# convex separable piecewise-linear penalty function.

###  PEOPLE DATA  ###

set Good;                              # Good risks
set Bad;                               # Bad risks
set people := Good union Bad;          # Everyone is either good or bad

param app_amt > 0;                     # General credit-approval amount
param bal_amt {people} >= app_amt;     # Maximum-balance amounts of individuals

###  FACTOR DATA  ###

set factors;                           # Individuals' responses to questions
set wt_types := {'pos','neg','free'};  # Required signs of weights
param wttyp {factors} symbolic in wt_types;
param answer {people,factors} >= 0;    # Numerical responses to all questions

###  PENALTY FUNCTION DATA  ###

# Parameters starting with G (or B) are for the good (or bad) risks

param Gpce > 1;                        # Linear pieces in penalty term
param Bpce > 1;
param Gslope {1..Gpce};  check {k in 1..Gpce-1}: Gslope[k] < Gslope[k+1];
param Bslope {1..Bpce};  check {k in 1..Bpce-1}: Bslope[k] < Bslope[k+1];
                                       # Increasing slopes in penalty term

set bkpt_types := {'A','B'};
param Gbktyp {1..Gpce-1} symbolic in bkpt_types;
param Bbktyp {1..Bpce-1} symbolic in bkpt_types;

param Gbkfac {1..Gpce-1};  check {k in 1..Gpce-2}: Gbkfac[k] <= Gbkfac[k+1];
param Bbkfac {1..Bpce-1};  check {k in 1..Bpce-2}: Bbkfac[k] <= Bbkfac[k+1];
                                       # Information to define the
                                       # increasing breakpoints in penalty
                                       # terms (see objective function)
param Gprop > 0;                       # Scale objective to simulate ratio
param Bprop > 0;                       # Gscale/Bscale of goods to bads
param Gratio := (Gprop / (Gprop+Bprop)) / 
                ((sum {i in Good} 1) / (sum {i in people} 1));
param Bratio := (Bprop / (Gprop+Bprop)) / 
                ((sum {i in Bad} 1) / (sum {i in people} 1));

###  VARIABLES  ###

var Wt_const;                          # Constant term in computing all scores
var Wt {j in factors} >= if wttyp[j] = 'pos' then 0 else -Infinity
                      <= if wttyp[j] = 'neg' then 0 else +Infinity;
                                       # Weights on the factors
var Sc {i in people};                  # Scores for the individuals

###  OBJECTIVE  ###

minimize penalty:                      # Sum of penalties for all individuals
   Gratio * sum {i in Good} << {k in 1..Gpce-1} if Gbktyp[k] = 'A' 
                                   then Gbkfac[k]*app_amt 
                                   else Gbkfac[k]*bal_amt[i];
                               {k in 1..Gpce} Gslope[k] >> Sc[i] +
   Bratio * sum {i in Bad}  << {k in 1..Bpce-1} if Bbktyp[k] = 'A' 
                                   then Bbkfac[k]*app_amt 
                                   else Bbkfac[k]*bal_amt[i];
                               {k in 1..Bpce} Bslope[k] >> Sc[i];

###  CONSTRAINTS  ###

def_Sc {i in people}:                  # Score = sum of answers times weights
   Sc[i] = Wt_const + sum {j in factors} answer[i,j] * Wt[j];
