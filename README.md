# RRmultinom Package

This package is designed to get the prevalence estimates of analyzing the sensitive randomized response questions `ever` and `last year` together. The package consists of three main functions.

## The function `prev_bivar`

It provides estimates of three categories: never, former and last year users, as well as a goodness-of-fit test statistic.

```{r}
 prev_bivar(Q1 = Q9, Q2 = Q10, data = Doping, p00 = 5/6, p11 = 5/6)
```

## The function \`rr_multinom\`

It fits a multinomial logistic regression model to the bivariate analysis of the "ever" and "last year" questions to investigate potential covariates associated with the sensitive characteristic.

```{r}
regression <- rr_multinom(formula = ~BodyBuild+Age, data = Doping,
                          response = c("Q9", "Q10"), p00 = 5/6, p11 = 5/6)
```

## The function \`AME\`

It calculates average marginal effect (AME) along with its standard error for all numerical and categorical covariates of the model object \\code{\\link{rr_multinom}}.

```{r}
ME <- AME(model = regression, data = Doping, method = "delta")
```

The AME can be calculated for subgroups not for the whole sample

```{r}
ME <- AME(model = regression, data = subset(Doping, BodyBuild == 0),
           method = "delta")
```
