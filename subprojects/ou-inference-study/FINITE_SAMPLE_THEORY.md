# Stationary OU Monte Carlo inference: finite-sample arguments

This addition specializes existing Monte Carlo rank-test reasoning to OU
inference and constructs a conservative, continuous-parameter cell envelope.
It does not establish priority for a new general testing principle. It uses
[Dufour's Monte Carlo test framework](https://jeanmariedufour.github.io/Dufour_1995_MCT_W.pdf)
(published in 2006, DOI 10.1016/j.jeconom.2005.06.007) and the standard
[Pinsker inequality](https://web.stanford.edu/class/stats311/lecture-notes.pdf).

## Assumptions and inferential target

Observe $n+1$ equally spaced values of a **stationary**, noise-free Gaussian OU
with $\theta>0$, unknown $\mu\in\mathbb R$, and unknown $\sigma>0$.
The first state has its stationary distribution and future Brownian increments
are independent of it. The sampling interval and record length are fixed
independently of the observed values. All statements below include randomness
in the initial state and simulation banks; they are not conditional on $X_0$.
Ideal calculations assume exact arithmetic and independent exact Gaussian draws.

The conditional profile likelihood is used as a **statistic**. Its calibration
comes from the unconditional stationary sampling distribution, even though the
statistic does not include an initial-density term. The previous fitted-model
bootstrap fixes observed $X_0$ and estimates nuisance parameters, so its
inferential construction differs.

## Proposition 1: nuisance-free stationary statistic

Let $V=\sigma^2/(2\theta)$ and $Z_k=(X_k-\mu)/\sqrt V$. Put
$\phi=e^{-\theta\Delta t}$. Then

$$Z_0\sim N(0,1),\qquad
Z_{k+1}=\phi Z_k+\sqrt{1-\phi^2}\,\xi_k,\qquad
\xi_k\stackrel{\mathrm{iid}}\sim N(0,1).$$

For the unrestricted regression, define $A$, $\widehat\phi$, $S_{\min}$, and
$S(\psi)=S_{\min}+A(\psi-\widehat\phi)^2$ as in [the profile derivation](THEORY.md).
Let $S_b=\inf_{0<\psi<1}S(\psi)$ and

$$T(X;\psi)=n\log\{S(\psi)/S_b\}.$$

**Proof.** Adding a constant to every state changes only the regression
intercept. Multiplying states by a nonzero constant leaves the slope unchanged
and multiplies both residual sums of squares by its square. Thus $T$ is affine
invariant, and $T(X;\psi)=T(Z;\psi)$. At a specified null slope, its stationary
law depends only on that slope and $n$, not on $\mu$ or $V$. The canonical
recursion above therefore simulates that null law without fitted nuisance
values. Nonconstant regressors and positive residual sums hold almost surely
for $n\ge3$ under this nonsingular Gaussian model. ∎

## Proposition 2: point-null rank test

For a specified $\psi\in(0,1)$, independently simulate $B$ canonical stationary
records under $\psi$. Evaluate $T_0=T(X;\psi)$ and $T_j$ on each simulated
record, and define

$$p_B=\frac{1+\sum_{j=1}^B\mathbf1\{T_j\ge T_0\}}{B+1}.$$

For any $0\le t\le1$,

$$P_\psi(p_B\le t)\le\frac{\lfloor(B+1)t\rfloor}{B+1}\le t.$$

**Proof.** Under the null, observed and simulated statistics are exchangeable.
With no ties, the descending rank of the observed statistic is uniform on
$\{1,\ldots,B+1\}$. Counting ties using $\ge$ gives a p-value at least as
large as an exchangeably tie-broken rank. This preserves the upper size bound.
The proof includes simulation randomness and needs no $n\to\infty$ argument. ∎

Consequently inverting tests with $p_B(\psi)>\alpha$ over the entire continuous
parameter space has coverage at least $1-\alpha$. Testing only grid points
does **not** supply that continuous confidence set. Piecewise interpolation
of p-values can exclude the true parameter between tested points.

## Proposition 3: exact stationary path KL divergence

Let $P_\phi$ denote the canonical law of $Z_0,\ldots,Z_n$, and let the reference
$\psi\in[0,1)$. The boundary reference zero is independent Gaussian white noise.
For $\phi\in[0,1)$,

$$D(P_\phi\Vert P_\psi)=\frac n2\left[
\log\frac{1-\psi^2}{1-\phi^2}
+\frac{1-\phi^2+(\phi-\psi)^2}{1-\psi^2}-1\right].\tag{1}$$

**Proof.** Initial densities are identical $N(0,1)$. Apply the KL chain rule
to the $n$ transition densities. The Gaussian conditional KL at preceding
state $z$ is half the log variance ratio plus half the variance and squared
mean-displacement ratio, minus one half. Under the stationary true model,
$E_\phi[Z_k^2]=1$, yielding (1) for each transition. ∎

For fixed $\psi$, differentiation gives

$$\frac{\partial D}{\partial\phi}
=n\left[\frac{\phi}{1-\phi^2}-\frac{\psi}{1-\psi^2}\right].$$

Since $x/(1-x^2)$ is strictly increasing on $[0,1)$, divergence decreases
toward its minimum at $\psi$ and increases after it. Its supremum on a closed
cell containing $\psi$ is therefore at an endpoint. If the cell reaches one,
the supremum diverges; its TV bound is safely set to one.

## Theorem: continuous-cell confidence envelope

Partition $[0,1)$ into a finite set of data-independent cells. For cell $I_j$,
choose a reference $\psi_j\in I_j$, and a bound

$$\delta_j\ge\sup_{\phi\in I_j}\|P_\phi-P_{\psi_j}\|_{TV}.$$

By Pinsker and Proposition 3, an available choice is

$$\delta_j=\min\left(1,\sqrt{\tfrac12\sup_{\phi\in I_j}
D(P_\phi\Vert P_{\psi_j})}\right).$$

If $\delta_j\ge\alpha$, retain the entire cell without testing. Otherwise run
the point-null rank test at $\psi_j$ at level $\alpha-\delta_j$, and retain
the entire cell if it is not rejected. Let $C_\phi$ be the union of retained
cells, and map it to $C_\theta$ by $\theta=-\log\phi/\Delta t$.
Then, for every true positive $\theta$, unknown $\mu$, and positive $\sigma$,

$$P_{\theta,\mu,\sigma}(\theta\in C_\theta)\ge1-\alpha.$$

**Proof.** Let $I_{j_0}$ contain the true slope. It is retained surely if its
bound is at least $\alpha$. Otherwise, under the reference distribution,
rejection probability is at most $\alpha-\delta_{j_0}$ by Proposition 2.
The simulation bank has the same law under true and reference data-generating
models and is independent of the observed record. Adding that bank to both
probability spaces does not increase their TV distance. Rejection under the
true model is therefore at most $\alpha-\delta_{j_0}+\delta_{j_0}=\alpha$.
Retaining the true cell includes the true parameter. Proposition 1 removes
the unknown location and scale. ∎

No Bonferroni factor is needed for coverage of the single true parameter:
only its deterministic cell enters the proof. Cells may use dependent banks,
but each bank must be independent of the observed record. For the empirical
coverage study, banks are independently regenerated for **every outer record**,
so binomial MC uncertainty calculations apply. A shared fixed bank could induce
dependence between outer outcomes and would require a different variance estimate.

## Implementation and limits

- The grid uses $u=\operatorname{atanh}\phi$, a fixed maximum $u=5$, and cell
  step no greater than $0.04/\sqrt n$. The first reference is zero; the last
  cell from $\tanh(5)$ to one is always retained. No small-theta tail is discarded.
- The returned set can be disconnected. Its enclosing interval is an explicitly
  conservative hull. The hull always reaches theta zero because of the retained
  tail; it cannot be advertised as a sharp two-sided interval.
- The coverage experiment tests membership in the true parameter's cell with
  independent banks. It does not compute full inversions or widths on every
  outer record. Full inversions are saved for predetermined replication-zero
  examples, with every cell decision retained.
- The theoretical KL bound is analytic. Code uses stable arithmetic, a small
  outward cushion, and independent covariance-matrix checks. This is not an
  interval-arithmetic computer certification for all floating-point inputs.
- Nonstationary starts, measurement errors, irregular sampling, or dependence
  between data and banks invalidate assumptions. Stress cases are diagnostics,
  not extensions of the theorem. High observed stress-case coverage is not a proof.
- Monte Carlo rank testing and Pinsker are established tools. Novelty of this
  particular OU cell construction has not been established by a systematic review.
