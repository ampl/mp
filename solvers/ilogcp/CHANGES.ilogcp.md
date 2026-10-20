Summary of recent updates to ilogcp for AMPL
============================================


## unreleased
- Options renamed to the MP style with the previous names kept as
  synonyms, e.g. *lim:time* (timelimit), *lim:sol* (solutionlimit),
  *lim:fail* (faillimit), *tech:threads* (workers), *tech:seed*
  (randomseed), *tech:outlev* (logverbosity, outlev), *mip:gap*
  (relativeoptimalitytolerance), *mip:gapabs* (optimalitytolerance),
  *tech:writemodel* (exportfile), *cp:searchtype* (searchtype), ...
- When the CPLEX optimizer is used, *lim:time*, *lim:sol*,
  *tech:threads*, *tech:seed*, *mip:gap*, *mip:gapabs* and *tech:outlev*
  are now applied to CPLEX; previously only CP Optimizer received them.
- Solve result codes follow the MP convention: limits return 400-series
  codes (time 402/472, solution limit 408, interrupt 401/471, fail limit
  430/490, branch limit 431/491, choice point limit 432/492),
  "infeasible or unbounded" returns 450; with CPLEX, the CPLEX status
  (time, node, iteration, memory, objective limits, numerical issues) is
  reported. See the `-!` output.
- New options: *tech:kpidisplay*, *tech:printmodeldetails*,
  *tech:modelanonymizer*, *tech:usefilelocations*, *tech:logsearchtags*,
  *cp:fds*, *cp:fdsmaxmemory*, *cp:searchtype=neighborhood*.
- IIS: *iis:find* runs the conflict refiner of CP Optimizer or CPLEX on
  infeasible problems and returns suffix .iis (mem/pmem/non) on
  constraints (and on variables with CP Optimizer); conflict refiner
  limits *lim:conflicttime*, *lim:conflictiter*, *lim:conflictbranch*,
  *lim:conflictfail* and *iis:onvariables*.
- Integer limit options are 64-bit (no more -1 shown for "no limit");
  *cp:temporalrelaxation* accepts ``auto`` (its default);
  corrected defaults in the option descriptions (time mode is elapsed
  time, optimality tolerance 1e-9, restart growth factor 1.15,
  log period 1000).


## 20260527
- *Linked with CPLEX and CP optimizer 22.2.0*


## 20250422
- *Linked with CPLEX and CP optimizer 22.01.02*


## 20230228
- *Linked with CPLEX and CP optimizer 22.01.01*


## 20221125
- *Linked with CPLEX and CP optimizer 22.01*
