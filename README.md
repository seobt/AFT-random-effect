# AFT-random-effect

1. Required data format
  - Data should be maxtrix
  - Firt column: Id number (e.g. patient id, center id)
  - Second column: log(survival time)
  - From third column: covariates
  - Last column: censoring indicator (1: observed, 0: censored)

2. function "SAFT" returns finite dimensional parameters and random effect distribution
3. function "PAFT" returns all parameters from random effect AFT model under normal random effect
4. For further information, please see "example.m"
   
