if (file.exists("_options.R")) source("_options.R")
library(panelPomp,quietly=TRUE)

TESTS_PASS <- NULL
test <- function(expr1,expr2,all="TESTS_PASS",env=parent.frame(),...)
  panelPomp:::test(expr1,expr2,all=all,env=env,...)

ppo <- panelRandomWalk(U=2,N=5)

pmf <- mif2(ppo,Np=10,rw.sd=rw_sd(X.0=0.2),cooling.fraction.50=0.5,
            cooling.type="geometric")

test(wQuotes(
  "Error : in ''traces'': name(s) ''par'' correspond to no parameter(s).\n"),
  traces(pmf,pars="par"))

## check whether all tests passed
all(get(eval(formals(test))$all))
if (!all(get(eval(formals(test))$all))) stop("Not all tests passed!")
