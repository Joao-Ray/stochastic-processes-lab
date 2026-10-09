# OU calibration: theoretical derivations

This note derives the formulas used in the [subproject](README.md). The
[numerical verification report](results/verification/report.md) evaluates them
against quadrature, independent optimization, recursive transitions, and Monte
Carlo experiments. A numerical check supports an implementation; the arguments
below justify the mathematical results under the stated assumptions.

Assume $\theta>0$, $\sigma>0$, constant $\mu$, and a standard Brownian motion.
Observations for estimation are equally spaced by $\Delta t>0$ and have no
measurement error. We condition on the initial observation rather than assume
its distribution is stationary.

## 1. Solve the SDE using an integrating factor

Let $Y_t=X_t-\mu$. Then

$$dY_t=-\theta Y_t\,dt+\sigma\,dW_t.$$

Since $e^{\theta t}$ is deterministic and has finite variation, the stochastic
product rule has no quadratic-covariation term involving it:

$$d(e^{\theta t}Y_t)=\theta e^{\theta t}Y_tdt+e^{\theta t}dY_t
=\sigma e^{\theta t}dW_t.$$

Integrate from $t$ to $t+h$ and multiply by $e^{-\theta(t+h)}$:

$$X_{t+h}=\mu+e^{-\theta h}(X_t-\mu)
+\sigma\int_t^{t+h} e^{-\theta(t+h-s)}dW_s.\tag{1}$$

The integrand is deterministic. Brownian increments after $t$ are independent
of the information available at $t$, so the integral is independent of $X_t$,
has mean zero, and is Gaussian.

## 2. Derive the transition variance with Itô isometry

For a deterministic square-integrable function $g$, Itô isometry gives
$E[(\int g(s)dW_s)^2]=\int g(s)^2ds$. Applying it to (1),

$$v(h)=\sigma^2\int_t^{t+h}e^{-2\theta(t+h-s)}ds
=\sigma^2\int_0^h e^{-2\theta u}du
=\frac{\sigma^2}{2\theta}(1-e^{-2\theta h}).\tag{2}$$

Therefore, conditional on $X_t=x$,

$$X_{t+h}\mid X_t=x\sim N(m(h;x),v(h)),\qquad
m(h;x)=\mu+e^{-\theta h}(x-\mu).\tag{3}$$

At equally spaced times this becomes

$$X_{k+1}=c+\phi X_k+\varepsilon_k,\qquad
\phi=e^{-\theta\Delta t},\quad c=\mu(1-\phi),\quad
q=\frac{\sigma^2}{2\theta}(1-\phi^2).\tag{4}$$

Innovations on disjoint intervals are independent $N(0,q)$. The observations
themselves are correlated; independent innovations do not imply independent
observations. `simulate_ou_exact` uses (4) on the grid.

## 3. Stationarity, autocorrelation, and half-life

Under stationarity let $V=\operatorname{Var}(X_k)$. From (4),

$$V=\phi^2V+q,\qquad V=\frac{q}{1-\phi^2}=\frac{\sigma^2}{2\theta}.$$

The invariant Gaussian law is $N(\mu,\sigma^2/(2\theta))$: substituting its mean
and variance into a Gaussian transition returns the same law.

Iterating (4), future innovations are independent of $X_k$, giving
$\operatorname{Cov}(X_{k+j},X_k)=\phi^jV$. Hence stationary autocorrelation at
elapsed lag $h$ is $e^{-\theta h}$. These statements require a stationary
initial law; the subproject's deterministic initial value is not stationary.

The expected displacement from $\mu$ halves when
$e^{-\theta h}=1/2$, so $h_{1/2}=\log(2)/\theta$.

## 4. Factor the conditional likelihood

For observed values $x_0,\ldots,x_n$, the Markov property factors the density
conditional on $x_0$ into $\prod_{k=0}^{n-1}p(x_{k+1}\mid x_k)$.
Taking logarithms of the Gaussian transitions gives

$$L(c,\phi,q)=-\ell
=\frac n2\log(2\pi q)+\frac{S(c,\phi)}{2q},\qquad
S(c,\phi)=\sum_{k=0}^{n-1}(x_{k+1}-c-\phi x_k)^2.\tag{5}$$

The conditioning is essential. Adding a stationary density for $x_0$ would
change this objective and generally change the estimator.

## 5. Derive the regression estimates and variance denominator

Put $u_k=x_k$, $y_k=x_{k+1}$, $r_k=y_k-c-\phi u_k$.
For fixed $q>0$, minimizing (5) is equivalent to minimizing $S$.
Setting its two derivatives to zero gives the normal equations

$$\sum r_k=0,\qquad\sum u_kr_k=0.$$

If $\sum(u_k-\bar u)^2>0$, solve them to obtain

$$\widehat\phi=\frac{\sum(u_k-\bar u)(y_k-\bar y)}{\sum(u_k-\bar u)^2},
\qquad\widehat c=\bar y-\widehat\phi\bar u.\tag{6}$$

The Hessian of $S$ is $2D^TD$ for design columns $(1,u)$, positive definite
when the regressors are nonconstant. Thus (6) is the unique global least-squares
minimum. For $S>0$, differentiation in $q$ gives

$$\frac{\partial L}{\partial q}=\frac n{2q}-\frac S{2q^2}=0,
\qquad\widehat q=\frac{S(\widehat c,\widehat\phi)}n.\tag{7}$$

At this optimum $\partial^2L/\partial q^2=n/(2\widehat q^2)>0$.
Substituting $q=S/n$ into (5) gives the profile objective
$L_{\mathrm{prof}}=\frac n2[\log(2\pi S/n)+1]$, strictly increasing in
$S>0$. Consequently the least-squares minimum also gives the global joint
minimum in $(c,\phi,q)$ when its residual sum of squares is positive.
The likelihood denominator is **$n$**, not $n-2$. Correlated lagged regressors
also mean the usual independent-design unbiasedness arguments cannot simply
be transferred to this estimator. We make no finite-sample unbiasedness claim.

## 6. Map the interior regression solution to OU parameters

If $0<\widehat\phi<1$ and $\widehat q>0$, invert (4):

$$\widehat\theta=-\frac{\log\widehat\phi}{\Delta t},\qquad
\widehat\mu=\frac{\widehat c}{1-\widehat\phi},\qquad
\widehat\sigma=\sqrt{\frac{2\widehat\theta\widehat q}{1-\widehat\phi^2}}.\tag{8}$$

This mapping is one-to-one between the interior OU parameter space and
$c\in\mathbb R$, $0<\phi<1$, $q>0$. Therefore the interior regression solution
also optimizes the conditional OU likelihood. If the unconstrained regression
slope is outside $(0,1)$, this argument does not produce an interior OU optimum.
`fit_ou` records that outcome as `OUFitError` rather than clipping the slope.

The mean estimate is especially sensitive when $\widehat\phi$ is near one,
because $\widehat\mu$ divides by $1-\widehat\phi$. High observation counts over
a short elapsed duration need not remove this difficulty.

## 7. Derive multi-step forecasts and composition

After $j$ steps,

$$X_{k+j}=\mu+\phi^j(X_k-\mu)
+\sum_{i=0}^{j-1}\phi^{j-1-i}\varepsilon_{k+i}.$$

Its conditional variance is $q\sum_{i=0}^{j-1}\phi^{2i}$, which equals (2)
at $h=j\Delta t$. More generally the transition moments compose as

$$m(h_1+h_2;x)=m(h_2;m(h_1;x)),
\quad v(h_1+h_2)=e^{-2\theta h_2}v(h_1)+v(h_2).\tag{9}$$

For known parameters and $h>0$, standardizing (3) gives a standard normal.
A central interval at level $a$ is $m\pm z_{(1+a)/2}\sqrt v$ and has exact
marginal coverage $a$. At $h=0$ the distribution is degenerate at $x$; coverage
of the zero-width interval is one. Plug-in parameters remove that exact-coverage
guarantee because estimation uncertainty is not included.

## 8. Derive the forecast-MSE benchmark

For a fixed prediction value $b$ at an origin $X_t=x$,

$$E[(X_{t+h}-b)^2\mid X_t=x]=v(h)+(m(h;x)-b)^2.\tag{10}$$

The cross term vanishes because the future noise has conditional mean zero.
For a predictor depending on the whole observed history, the same identity
holds conditional on that history, with $x=X_t$ and $b$ then known.
The conditional mean is therefore the squared-error optimal predictor, with
MSE $v(h)$. Persistence uses $b=x$ and has MSE
$v(h)+(\mu-x)^2(1-e^{-\theta h})^2$. This conditional benchmark fixes $x$;
it is distinct from averaging over random forecast origins in `run.py`.

For Monte Carlo error bars, if $Y\sim N(d,v)$ then
$\operatorname{Var}(Y^2)=2v^2+4d^2v$. Consequently an average of $M$ independent
squared errors has standard error $\sqrt{(2v^2+4d^2v)/M}$.

## 9. Derive the paired mean-shift and lost-coverage formulas

Suppose two futures start at the same value and share Brownian noise, but one
equilibrium mean changes by $\delta$. Subtract their SDEs. Their difference
$D_h$ satisfies the deterministic ODE

$$D'_h=\theta(\delta-D_h),\qquad D_0=0,\qquad
D_h=\delta(1-e^{-\theta h}).\tag{11}$$

This proves the paired-shift construction in `run.py`. The shifted variance
is still $v(h)$. Relative to the old known-parameter forecast mean, the shifted
future has displacement $D_h$. Write $r=z_{(1+a)/2}\sqrt{v(h)}$. Its coverage
under the old interval is

$$\Phi\!\left(\frac{r-D_h}{\sqrt{v(h)}}\right)
-\Phi\!\left(\frac{-r-D_h}{\sqrt{v(h)}}\right).\tag{12}$$

An oracle using the correct shifted mean restores coverage $a$. Equation (12)
assumes known old parameters; it is not the exact coverage formula for a fitted
plug-in interval. We verify that simpler benchmark separately.

## 10. Theory-to-verification map

| Derived statement | Independent numerical comparison | Recorded quantity |
| --- | --- | --- |
| Itô variance, (2) | Adaptive quadrature of the squared kernel | Absolute and relative error |
| Conditional moments, coverage and MSE, (3), (10) | 50,000 simulated paths from a fixed initial value | Moment errors, Monte Carlo SEs, coverage counts |
| Interior conditional MLE, (5)–(8) | Two-start numerical optimization in physical parameters on eight seeded records | Objective gap, parameter discrepancy, convergence status |
| Multi-step composition, (9) | Recursive mean/variance transitions | Maximum discrepancy |
| Paired shift, (11) | Coupled exact recursions plus numerical ODE solution | Path-difference and ODE errors |
| Old-interval shifted coverage, (12) | Simulated shifted futures | Empirical and analytic coverage, Monte Carlo SE |

Stationarity and autocorrelation are derived here and investigated in the
parent OU notebook; this new verification focuses on calibration and prediction.
Monte Carlo agreement is not a proof. Scientific 95% intervals are reported
without selecting favorable seeds; automated six-SE checks detect gross
distributional implementation errors and are not 95% significance tests.

Background: [Oxford's stochastic calculus course](https://courses.maths.ox.ac.uk/course/info.php?id=968)
covers stochastic integration and Itô isometry;
[Berkeley's time-series lecture on conditional AR likelihood](https://www.stat.berkeley.edu/~aditya/resources/LectureTHIRTEEN.pdf)
explains the least-squares connection. The model-specific steps above are
derived directly for the SDE and conditioning used in this repository.
