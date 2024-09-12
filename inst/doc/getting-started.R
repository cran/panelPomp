## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(panelPomp)

## ----createPpomp1-------------------------------------------------------------
mod1 <- pomp::gompertz()  # Using default values 
mod2 <- pomp::gompertz(K = 2, r = 0.01)  # Overwriting some of the defaults 
mod3 <- pomp::gompertz(K = 1.5, sigma = 0.15, X_0 = 5)
mod4 <- pomp::gompertz(K = 1.5, r = 0.05, X_0 = 5)
mod5 <- pomp::gompertz(K = 5, sigma = 0.08)

panelMod1 <- panelPomp(
  object = list(mod1, mod2, mod3, mod4, mod5)
)

## ----showMod1-----------------------------------------------------------------
print(specific(panelMod1))
print(shared(panelMod1))

## ----createPpomp2-------------------------------------------------------------
panelMod2 <- panelPomp(
  object = list(mod1, mod2, mod3, mod4, mod5),
  shared = c("tau" = 0.1)
)

specific(panelMod2)
shared(panelMod2)

## ----unitSpecific-vector------------------------------------------------------
specific(panelMod2, format = 'vector')

## ----setterFuns---------------------------------------------------------------
shared(panelMod2) <- c('r' = 0.05, 'sigma' = 0.1)
specific(panelMod2) <- c('tau[unit1]' = 0.11, 'tau[unit4]' = 0.09)

print(shared(panelMod2))
print(specific(panelMod2))

## -----------------------------------------------------------------------------
specific(panelMod2) <- matrix(
  data = rbind(c(1.24, 1.78), 
               c(   2,    3)),
  nrow = 2, 
  dimnames = list(
    param = c("K", "X_0"), 
    unit = c('unit2', 'unit4')
  )
)

specific(panelMod2)

## ----errorShared, error=TRUE--------------------------------------------------
shared(panelMod2) <- c("foo" = 1)

## ----errorSpecific, error=TRUE------------------------------------------------
specific(panelMod2) <- c("tau[unit6]" = 1)

## ----rinit--------------------------------------------------------------------
rinit_u <- pomp::Csnippet(
"
X = X_0;
"
)

## ----rprocess-----------------------------------------------------------------
rprocess_u <- pomp::Csnippet(
"
double S = exp(-r);
double eps = exp(rnorm(0,sigma));

X = pow(K, 1-S) * pow(X, S) * eps;
"
)

## ----rmeasure-----------------------------------------------------------------
rmeasure_u <- pomp::Csnippet(
"
Y = rlnorm(log(X),tau);
"
)

## ----dmeasure-----------------------------------------------------------------
rmeasure_u <- pomp::Csnippet(
"
lik = dlnorm(Y,log(X),tau,give_log);
"
)

