# Derivations for the OU uncertainty benchmark

The parent [OU derivation](../ou-calibration/THEORY.md) establishes the exact
transition and conditional MLE. Here mu and sigma remain unknown.

## 1. Dimensionless experiment design

Let V=sigma²/(2 theta), s=theta*t, and Z_s=(X_(s/theta)-mu)/sqrt(V).
Brownian scaling gives dZ_s=-Z_s ds+sqrt(2)dB_s. With a stationary initial law,
Z_0~N(0,1). Thus normalized experiments are indexed by r=theta*T and
a=theta*dt, with n=r/a transitions. More observations at fixed r do not
increase the number of elapsed reversion times. A scale-invariance test
checks the implementation under a change of time units.

## 2. Exact conditional profile likelihood

Write u_k=x_k, y_k=x_(k+1), and center both vectors. For a fixed slope phi,
the nuisance intercept optimum is c(phi)=mean(y)-phi*mean(u). Define

$$A=\sum(u_k-\bar u)^2,\quad
\widehat\phi=\frac{\sum(u_k-\bar u)(y_k-\bar y)}A.$$

Expanding the residual square and using the normal equations gives

$$S(\phi)=S_{\min}+A(\phi-\widehat\phi)^2.$$

Profiling the innovation variance q gives q(phi)=S(phi)/n. Because
q=sigma²(1-phi²)/(2 theta), optimizing freely over q at a fixed interior phi
is equivalent to optimizing sigma>0. Likewise c is free because mu=c/(1-phi).

On the physical interval 0<phi<1 the infimum of S is at the unconstrained
slope when interior, otherwise at a boundary limit. Put
phi_b=min(max(phi_hat,0),1) solely for calculating this infimum. It is not a
reported OU point estimate: phi=1 or 0 corresponds to a boundary supremum,
possibly with divergent nuisance parameters. Define S_b=S(phi_b).

The profile likelihood-ratio statistic is

$$\Lambda(\phi)=n\log\{S(\phi)/S_b\}.$$

For a cutoff c_alpha>=0, Lambda<=c_alpha is equivalent to

$$|\phi-\widehat\phi|\le
\sqrt{\{S_b e^{c_\alpha/n}-S_{\min}\}/A}.$$

Intersect these limits with (0,1), then transform by theta=-log(phi)/dt.
An upper phi limit of 1 yields a theta lower limit of zero; a lower phi
limit of zero yields an infinite theta upper limit. These are closures for
reporting; positive truth values have unchanged membership. Boundary cases
are not forcibly turned into successful point fits.

The chi-square(1) cutoff is a **regular asymptotic approximation**, not an exact
small-sample guarantee. Boundary suprema and weak identification can invalidate
the regular Wilks approximation; their coverage is measured, not assumed.

## 3. Wald interval and the sampling-density distinction

The conditional regression Hessian gives Var(phi_hat) approximately q/A.
The delta derivative dtheta/dphi=-1/(dt*phi) gives

$$\widehat{\mathrm{SE}}(\widehat\theta)
=\frac{\sqrt{\widehat q/A}}{\Delta t\widehat\phi}.$$

The interval is theta_hat +/- z_.975*SE, intersected with theta>=0.
This approximation is not an exact regression t interval with lagged responses.

For stationary expanding-domain OU observations, replacing A by nV leads to
the heuristic variance (exp(2a)-1)/(n dt²). At a->0 and n=T/dt it approaches
2theta/T rather than zero when T is fixed. The substitution is an
**expanding-domain heuristic**, not a fixed-T asymptotic theorem: realized A
remains random at fixed T. This distinction motivates the paired infill study.

## 4. Bootstrap procedures and correction

Given an interior fit, simulate B exact fitted-OU training paths of identical
length with initial state equal to the observed x_0. Refit each in unrestricted
AR coordinates. Retain all phi_hat_star values, including those outside (0,1).

If Q_p denotes the empirical slope quantile, the basic bootstrap interval is
[2phi_hat-Q_.975, 2phi_hat-Q_.025]. Intersect with (0,1) and transform. An
empty intersection is an empty confidence set, not missing information.

For bootstrap profile calibration, calculate Lambda_star(phi_hat) on each
bootstrap record, where its denominator uses that record's constrained
likelihood supremum. The empirical 95th percentile replaces the chi-square
cutoff in Section 2. The fitted nuisance parameters remain plug-in values;
this is not exact finite-sample pivotal inversion.

The optional slope correction phi_bc=2phi_hat-mean(phi_hat_star) estimates
and subtracts bootstrap slope bias. Map to theta only when 0<phi_bc<1.
By nonlinearity, this is not equivalent to correcting theta directly, nor
does reducing bias guarantee smaller RMSE. The main comparisons use common
admissible records and report correction failures separately.

## 5. Why observation error defeats a noise-free bootstrap

For stationary latent OU with variance V, observe Y_k=X_k+e_k with independent
e_k~N(0,eta²). Then Var(Y_k)=V+eta² and lag-one covariance is phi*V.
At fixed sampling interval and expanding duration the regression pseudo-slope is

$$\phi_{\mathrm{pseudo}}=\frac{\phi V}{V+\eta^2},\qquad
\theta_{\mathrm{pseudo}}=\theta+\frac{\log(1+\eta^2/V)}{\Delta t}.$$

This is an ergodic long-record limit, not the exact expectation of a finite
sample estimate. The lag sequence of Y is generally not that of an AR(1),
so independent OU innovations are misspecified. Resampling a fitted noise-free
OU cannot correct this structural mismatch. The stress experiment targets
latent theta; pseudo-theta is reported separately, not substituted as truth.

## 6. Monte Carlo uncertainty and comparisons

For M independent outer attempts, coverage p_hat has estimated MC SE
sqrt(p_hat(1-p_hat)/M); Wilson intervals accompany small/near-boundary counts.
Conditional coverage uses the number of available sets instead of M.
An unavailable method contributes zero only to the explicitly named
success-and-coverage measure.

For errors e_i, bias SE is sd(e_i)/sqrt(M). RMSE=sqrt(mean(e_i²)); its
delta-method MC SE is sd(e_i²)/(2*RMSE*sqrt(M)). Both are conditional on
the stated admissibility denominator. For paired method hits, compute
sd(hit_A-hit_B)/sqrt(M_common), not a sum of independent coverage variances.
Intervals across sampling grids and paired noise levels are dependent.

Finite mean widths condition on finite nonempty sets and are accompanied by
their counts. If any width is infinite, the unconditional sample mean width
is infinite; it is never represented as a finite mean. No multiplicity-
adjusted superiority claims are made from the exploratory comparisons.
