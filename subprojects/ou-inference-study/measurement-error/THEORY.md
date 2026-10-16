# Noisy stationary OU: likelihood, replication and information limits

These are explicit model-specific derivations using established Gaussian
likelihood, REML, Fisher information and two-point testing tools. Priority or
novelty is not claimed. Calibration of profile likelihood below is asymptotic.

## 1. Observation law and identifiability

At site i with k simultaneous readings, let

$$Y_{ij}=X_{i\Delta}+\epsilon_{ij},\quad
\epsilon_{ij}\stackrel{iid}{\sim}N(0,\eta^2),\quad
X\text{ stationary OU},\quad V=\sigma^2/(2\theta).$$

All measurement errors are independent of X. Put phi=exp(-theta*Delta),
s^2=V+eta^2, gamma=V/s^2. The raw covariance is

$$s^2[\gamma R_\phi\otimes J_k+(1-\gamma)I_{mk}],\qquad
(R_\phi)_{ij}=\phi^{|i-j|}.$$

For one reading per site, population covariance satisfies C0=V+eta^2,
C1=V*phi, C2=V*phi^2. If V>0 and 0<phi<1, these identify
phi=C2/C1, V=C1^2/C2, eta^2=C0-V. This population statement does not make
the sample ratios stable, admissible or unbiased. At gamma=0 there is no
latent signal and the observation law is independent of drift.

Noisy Y is not a first-order Markov process. Subtracting phi*Y_{i-1} gives
latent innovation + epsilon_i - phi*epsilon_{i-1}; adjacent residuals have
covariance -phi*eta^2. A product of one-step noisy pair densities generally
is not the full joint likelihood. We implement the joint covariance likelihood,
checked against sequential filtering, rather than treating observations as AR(1).

## 2. Exact reduction for repeated readings

Let barY_i be site means, W=sum_ij(Y_ij-barY_i)^2 and d=m(k-1).
An orthogonal transformation within each site separates the scaled mean
sqrt(k)*barY_i and k-1 independent contrasts. Gaussian orthogonality implies

$$\bar Y\sim N(\mu\mathbf1,s^2 A),\quad
A=\gamma R_\phi+\frac{1-\gamma}{k}I_m,\qquad
\frac{W}{\eta^2}\sim\chi^2_d,$$

and W is independent of all site means. Thus within-site contrasts contain
information about sensor variance without containing latent fluctuations.
Replicates do not create additional independent latent states.

The raw shape log determinant is

$$L=\log|A|+m\log k+d\log(1-\gamma).$$

For k=1 omit the last term. GLS profiles the mean as

$$\widehat\mu=\frac{\mathbf1^TA^{-1}\bar y}{\mathbf1^TA^{-1}\mathbf1},\quad
S=(\bar y-\widehat\mu\mathbf1)^TA^{-1}(\bar y-\widehat\mu\mathbf1)
+\frac{W}{1-\gamma}.$$

For ML put nu=N=mk. For REML put nu=N-1. The profiled objective, up to fixed
constants, is

$$\ell_{ML}=\tfrac12[L+N\{\log(2\pi S/N)+1\}],$$

$$\ell_{REML}=\tfrac12[L+\log(\mathbf1^TA^{-1}\mathbf1)
+(N-1)\{\log(2\pi S/(N-1))+1\}].$$

**Derivation.** The mean GLS equation minimizes the quadratic; the residual
contrasts supply W/(1-gamma). Optimizing scale gives S/nu. For REML integrate
or remove the one-dimensional fixed mean, contributing the log mean-information
term and removing one residual degree of freedom. Orthonormal residual
contrasts give the same objective minus the constant log(N)/2. Scale therefore
is not generally an unbiased finite-sample estimator after optimizing shape.

For replicated data with positive W, gamma=1 is singular and has zero likelihood.
With W=0 the boundary is a degenerate Gaussian measure, requiring a separate
zero-noise analysis; the public fitter rejects this singular replicated input.
The replicated experiments all have strictly positive sensor variance.

## 3. Linear-time covariance evaluation

Write q=1-phi^2 and Q=q*R_phi^{-1}. Q is tridiagonal: endpoints 1,
interior diagonal 1+phi^2, and off diagonal -phi. Define

$$H=\gamma q I+\frac{1-\gamma}{k}Q.$$

Then A=R_phi*H/q and, because Q and H commute,

$$A^{-1}b=H^{-1}Qb,\qquad \log|A|=\log|H|-\log q.$$

**Proof.** Multiply A by Q using R_phi*Q=qI. The AR correlation determinant
is |R_phi|=q^(m-1), so |A|=|R_phi|*|H|/q^m=|H|/q. Banded Cholesky of H
and two triangular solves evaluate the log determinant and GLS products in
O(m) work. Dense covariance, Kalman updates and residual-contrast references
provide independent implementation checks.

This algebra speeds an established likelihood; it is not claimed as a new
inference principle. Numerical multi-start optimization is not a global proof.
The physical parameter space is open; finite search caps and white-noise
limits are explicitly flagged rather than reported as successful finite drifts.

## 4. Local information at a fixed scalar-reading budget

For the averaged data, let C=V*R_phi+eta^2/k*I. With parameters
beta=(theta,V,eta^2), derivatives are

$$C_\theta=-V D\circ R_\phi,\quad C_V=R_\phi,\quad C_{\eta^2}=I/k,$$

where D_ij=Delta*|i-j|. Gaussian covariance-score information is

$$I_{ab}=\tfrac12\operatorname{tr}(C^{-1}C_a C^{-1}C_b)
+\mathbf1_{a=b=\eta^2}\frac{d}{2\eta^4}.$$

**Proof.** The expected product of centered Gaussian quadratic scores is the
trace term. Independent within-site contrasts have d degrees of freedom at
variance eta^2; their variance score supplies d/(2*eta^4). Mean scores are
orthogonal in expectation to covariance scores. Eliminating covariance
nuisances V and eta^2 gives the Schur complement

$$I_{eff}=I_{\theta\theta}-I_{\theta\lambda}I_{\lambda\lambda}^{-1}
I_{\lambda\theta},\quad \lambda=(V,\eta^2).$$

The local Cramer-Rao variance benchmark is 1/I_eff for regular unbiased
estimators at interior parameters. It does not lower-bound the RMSE of every
biased estimator, give small-sample CI coverage or ensure k=2 is always better.
At fixed N=m*k and T, increasing k trades time resolution for independent
sensor-variance information. Actual k comparisons use independent trajectories.

## 5. A finite-sample weak-signal limitation

Consider two slopes phi0,phi1, the same mean and unit total variance, and
signal fraction gamma. Their covariance matrices for m single readings are
P=gamma*R0+(1-gamma)I and Q=gamma*R1+(1-gamma)I. They are admissible OU
submodels for gamma>0; diffusion changes with drift to keep stationary V=gamma.
This restriction is enough for a lower bound over a larger unknown-variance class.

The exact Gaussian divergence is

$$D(P\Vert Q)=\tfrac12[\operatorname{tr}(Q^{-1}P)-m+
\log|Q|-\log|P|].$$

For 0<=gamma<1,

$$D(P\Vert Q)\le
\frac{\gamma^2\|R_0-R_1\|_F^2}{4(1-\gamma)^2},\qquad
\|R_0-R_1\|_F^2=2\sum_{h=1}^{m-1}(m-h)(\phi_0^h-\phi_1^h)^2.$$

**Proof.** Put f(M)=-log|M|. Gaussian KL is one half its Bregman divergence
f(P)-f(Q)-<grad f(Q),P-Q>. The Hessian on a direction E is
tr(M^{-1}E M^{-1}E)<=||E||_F^2/lambda_min(M)^2.
All matrices on the line from P to Q have minimum eigenvalue at least 1-gamma.
Taylor's integral remainder contributes a factor 1/2; the Gaussian prefactor
contributes another 1/2. Insert E=gamma*(R0-R1). The diagonal of R0-R1 is
zero; counting its off-diagonals gives the displayed sum. This bound is valid
for a fixed finite m; asymptotic rates require a specified m/gamma/sampling regime.

Let delta=min(1,sqrt(min(D(P||Q),D(Q||P))/2)). Pinsker implies TV(P,Q)<=delta.
Then:

- Any test with size <=alpha at phi0 has power at phi1 at most alpha+delta
  (capped at one), by the event definition of total variation.
- If a random confidence set covers each of the two positive drift values with
  probability >=1-alpha under its own law, then under P the probability of
  containing **both** is >=1-2*alpha-delta. Transfer coverage of the second
  value from Q to P using TV, then apply a union bound.
- For any drift estimator and separation Delta_theta, nearest-alternative
  classification gives max_i P_i(|theta_hat-theta_i|>=Delta_theta/2)
  >=(1-delta)/2. Consequently the maximum squared-error risk is at least
  Delta_theta^2*(1-delta)/8.

Independent simulation randomization can be appended to both laws without
increasing TV. These limits therefore also apply to randomized inference.
As gamma tends to zero, information about drift disappears at fixed m; no
uniformly narrow, powerful method can be rescued solely by changing optimizers.
This does not prove that any particular fitted CI has correct coverage, nor
that replication always improves drift estimation. The analytic bound and
exact covariance KL are separate from observed likelihood performance.

## 6. Relation to prior work and limits

[Chen, Simpson & Ying (2000)](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n17.pdf)
already study OU plus white noise, infill likelihood asymptotics and the
microergodic product V*theta. Their results must not be described as drift
consistency with both V and theta unknown on a fixed time window.
[Holý & Tomanová, arXiv:1811.09312](https://arxiv.org/abs/1811.09312)
already treat noisy OU estimation and ARMA representations.
[Calderon (2013)](https://arxiv.org/abs/1304.4196) studies finite-sample bias
and noise. REML is an established variance-component method. The information
and two-point arguments above are applications of standard Gaussian and
information-theoretic tools, not established original theorems.

Noisy profile-LR sets use an asymptotic chi-square cutoff; boundary mixtures,
short windows and weak signal can invalidate its nominal 95% interpretation.
The preceding noise-free rank theorem is not imported into this noisy model.
Only true-point membership is measured here. Full-set widths, uniform finite-
sample noisy coverage, irregular sampling, and colored-noise robustness remain
unresolved. The colored-noise/fixed-start runs are diagnostics, not extensions
of the assumed likelihood.
