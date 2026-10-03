# Mathematical Notes and Review Guide

These notes collect the main arguments behind the five numerical experiments. They are meant to be read together with the notebooks: the notes explain why the formulas are true, while the notebooks test what those formulas look like in finite simulations.

## How to read the experiments

Every experiment separates three kinds of statements:

1. **Model assumptions** define the random process, such as independent increments or the Markov property.
2. **Theoretical consequences** follow from those assumptions, such as a mean, variance, stationary distribution, or limiting law.
3. **Numerical evidence** checks whether a finite implementation behaves consistently with the theory.

A simulation can reveal a coding error or show agreement at a chosen scale. It cannot prove a theorem or establish that a real dataset follows the model.

## 1. Symmetric random walk

Let

$$
S_n=\sum_{i=1}^n X_i,
\qquad
P(X_i=1)=P(X_i=-1)=\frac12,
$$

where the increments are independent.

### Mean and variance

For one increment,

$$
E[X_i]=\frac12-\frac12=0,
\qquad
E[X_i^2]=1.
$$

Therefore

$$
\operatorname{Var}(X_i)=E[X_i^2]-E[X_i]^2=1.
$$

Linearity of expectation gives

$$
E[S_n]=\sum_{i=1}^n E[X_i]=0.
$$

Independence makes all covariance terms vanish, so

$$
\operatorname{Var}(S_n)
=\sum_{i=1}^n\operatorname{Var}(X_i)
=n.
$$

The standard deviation is therefore $\sqrt n$. A typical displacement grows like the square root of elapsed time even though the path contains $n$ unit-sized steps.

### Exact endpoint distribution

Let $B_n$ be the number of $+1$ steps. Then $B_n\sim\operatorname{Binomial}(n,1/2)$ and

$$
S_n=B_n-(n-B_n)=2B_n-n.
$$

Thus $S_n$ has the same parity as $n$, and for an admissible endpoint $k$,

$$
P(S_n=k)
=\binom{n}{(n+k)/2}2^{-n}.
$$

This formula explains the lattice of possible endpoints used in the notebook's distribution comparison.

### First hitting time

For a positive level $a$, define

$$
T_a=\inf\{n\geq0:S_n=a\}.
$$

The one-dimensional symmetric walk reaches $a$ eventually with probability one, but $E[T_a]=\infty$. These facts can coexist because the hitting-time distribution has a heavy tail. A finite simulation horizon therefore creates right-censored paths: an unfinished path has an unknown later hitting time, rather than no hitting time.

Here is a finite-interval argument for both claims. For any positive integer $b$, let

$$
\tau_b=\inf\{n\geq0:S_n\in\{a,-b\}\}.
$$

The martingale $S_n$ and absorption at either boundary give the gambler's-ruin probabilities

$$
P(S_{\tau_b}=a)=\frac{b}{a+b},
\qquad
P(S_{\tau_b}=-b)=\frac{a}{a+b}.
$$

Because hitting $a$ before $-b$ guarantees eventual hitting of $a$,

$$
P(T_a<\infty)\geq\frac{b}{a+b}.
$$

Letting $b\to\infty$ proves $P(T_a<\infty)=1$. The martingale $S_n^2-n$ gives

$$
E[\tau_b]=E[S_{\tau_b}^2]
=a^2\frac{b}{a+b}+b^2\frac{a}{a+b}
=ab.
$$

Since $\tau_b\leq T_a$, we have $E[T_a]\geq ab$ for every $b$. The bound grows without limit, so $E[T_a]=\infty$. The use of optional stopping here is justified for the walk absorbed between two finite boundaries; it does not assume a finite mean for $T_a$.

### Code connection

- `simulate_random_walk` constructs cumulative sums of independent $\pm1$ increments.
- `endpoint_statistics` computes the empirical mean and unbiased sample variance.
- `simulate_first_hitting_times` records the first visit and preserves censored observations. The notebook evaluates the exact endpoint probabilities with the binomial distribution from SciPy.

## 2. Homogeneous Poisson process

A rate-$\lambda$ Poisson process $N(t)$ starts from zero and has stationary, independent increments. For $0\leq s<t$,

$$
N(t)-N(s)\sim\operatorname{Poisson}(\lambda(t-s)).
$$

### Count moments

If $Z\sim\operatorname{Poisson}(m)$, then

$$
E[Z]=m,
\qquad
\operatorname{Var}(Z)=m.
$$

Taking $m=\lambda t$ gives

$$
E[N(t)]=\operatorname{Var}(N(t))=\lambda t.
$$

The equality of mean and variance is a characteristic consequence of the Poisson model, although matching these two moments alone does not identify a distribution.

### Waiting-time distribution

Let $W_1$ be the first arrival time. The event $W_1>t$ is exactly the event of no arrivals by time $t$. Hence

$$
P(W_1>t)=P(N(t)=0)=e^{-\lambda t},
$$

which is the survival function of $\operatorname{Exp}(\lambda)$. Therefore

$$
E[W_1]=\frac1\lambda,
\qquad
\operatorname{Var}(W_1)=\frac1{\lambda^2}.
$$

Independent, stationary increments imply that every later inter-arrival time has the same exponential law and is independent of the previous waits.

### Memorylessness

For $s,u\geq0$,

$$
P(W>s+u\mid W>s)
=\frac{e^{-\lambda(s+u)}}{e^{-\lambda s}}
=e^{-\lambda u}
=P(W>u).
$$

This is a conditional distribution statement. It does not say that an observed waiting time has no elapsed duration.

### Count and arrival-time duality

If $A_k=W_1+\cdots+W_k$ is the time of the $k$th arrival, then

$$
N(t)\geq k \quad\Longleftrightarrow\quad A_k\leq t.
$$

The count process and the sequence of waiting times are two descriptions of the same arrivals.

### Code connection

- `simulate_count_paths` samples independent Poisson increments on a time grid and cumulatively sums them; the notebook also uses those increments on non-overlapping intervals.
- `simulate_interarrival_times` samples exponential gaps.
- `simulate_arrival_times` cumulatively adds exponential gaps to construct jump times up to a fixed horizon.

## 3. Finite-state Markov chain

For a time-homogeneous Markov chain with transition matrix $P$,

$$
P_{ij}=P(X_{n+1}=j\mid X_n=i).
$$

Each entry is nonnegative and each row sums to one.

### Why matrix powers give multi-step probabilities

Condition on the intermediate state $k$:

$$
\begin{aligned}
P(X_2=j\mid X_0=i)
&=\sum_k P(X_1=k\mid X_0=i)
          P(X_2=j\mid X_1=k)\\
&=\sum_k P_{ik}P_{kj}\\
&=(P^2)_{ij}.
\end{aligned}
$$

Repeating the argument proves the Chapman--Kolmogorov relation

$$
P(X_n=j\mid X_0=i)=(P^n)_{ij}.
$$

If the initial distribution is the row vector $\alpha_0$, then $\alpha_n=\alpha_0P^n$.

### Stationary distribution

A stationary distribution $\pi$ satisfies

$$
\pi P=\pi,
\qquad
\sum_i\pi_i=1.
$$

For

$$
P=\begin{pmatrix}0.8&0.2\\0.4&0.6\end{pmatrix},
$$

the balance equations give

$$
0.2\pi_{\mathrm{Sunny}}=0.4\pi_{\mathrm{Rainy}},
$$

so

$$
\pi=\left(\frac23,\frac13\right).
$$

Stationarity describes stable population proportions. Individual paths continue moving between states.

### Two kinds of convergence

For this irreducible, aperiodic finite chain:

1. the distribution $\alpha_0P^n$ converges to $\pi$;
2. the fraction of time that one long path spends in each state converges to $\pi$.

The notebook displays both statements because they use different averages.

### Code connection

- `n_step_transition` computes $P^n$.
- `stationary_distribution` solves the balance equations plus normalization.
- `simulate_chains` samples the next state from the row associated with the current state.
- `empirical_transition_matrix` reconstructs transition probabilities from observed moves.

## 4. Ornstein--Uhlenbeck process

The OU process solves

$$
dX_t=\theta(\mu-X_t)\,dt+\sigma\,dW_t,
\qquad \theta>0.
$$

The drift points toward $\mu$, while the Brownian term continuously introduces random variation.

### Exact representation

Multiply $X_t-\mu$ by the integrating factor $e^{\theta t}$. Itô's product rule yields

$$
d\!\left[e^{\theta t}(X_t-\mu)\right]
=\sigma e^{\theta t}\,dW_t.
$$

Integrating from zero to $t$ gives

$$
X_t=\mu+(x_0-\mu)e^{-\theta t}
+\sigma\int_0^t e^{-\theta(t-s)}\,dW_s.
$$

The stochastic integral is Gaussian with mean zero. Itô isometry gives its variance:

$$
\sigma^2\int_0^t e^{-2\theta(t-s)}\,ds
=\frac{\sigma^2}{2\theta}\left(1-e^{-2\theta t}\right).
$$

Consequently,

$$
E[X_t]=\mu+(x_0-\mu)e^{-\theta t},
$$

$$
\operatorname{Var}(X_t)
=\frac{\sigma^2}{2\theta}\left(1-e^{-2\theta t}\right).
$$

### Stationarity, autocorrelation, and half-life

As $t\to\infty$,

$$
X_t\xrightarrow{d}N\!\left(\mu,\frac{\sigma^2}{2\theta}\right).
$$

In stationarity, the autocorrelation at lag $\tau$ is

$$
\rho(\tau)=e^{-\theta\tau}.
$$

The displacement of the conditional mean from $\mu$ is multiplied by $e^{-\theta t}$. Setting this factor to $1/2$ gives

$$
t_{1/2}=\frac{\log 2}{\theta}.
$$

Thus $\mu$ sets the equilibrium level, $\sigma$ controls shock size, and $\theta$ controls the reversion speed. With $\sigma$ fixed, increasing $\theta$ also reduces stationary variance.

### Euler--Maruyama approximation

For a time step $\Delta t$,

$$
X_{k+1}=X_k+\theta(\mu-X_k)\Delta t
+\sigma\sqrt{\Delta t}\,Z_k,
\qquad Z_k\sim N(0,1).
$$

The discrete recursion has autoregressive coefficient $1-\theta\Delta t$. The implementation requires $\theta\Delta t<2$ so that the deterministic recursion does not become unstable. A smaller step reduces discretization error but requires more computation.

### Code connection

- `ou_mean`, `ou_variance`, and `stationary_variance` evaluate exact formulas.
- `ou_autocorrelation` and `ou_half_life` quantify memory and reversion speed.
- `simulate_ou_euler` implements the Euler--Maruyama recursion.
- `sample_autocorrelation` estimates temporal dependence from one path.

## 5. Monte Carlo methods

Monte Carlo estimation rewrites a target as an expectation and approximates that expectation with a sample average.

### Estimating pi

For independent $(U_i,V_i)$ uniformly distributed on $[0,1]^2$, define

$$
Y_i=4\mathbf{1}\{U_i^2+V_i^2\leq1\}.
$$

The event inside the indicator is a quarter circle of area $\pi/4$. Therefore

$$
E[Y_i]=4\frac\pi4=\pi.
$$

Writing $p=\pi/4$ and using the variance of a Bernoulli variable,

$$
\operatorname{Var}(Y_i)
=16p(1-p)
=4\pi-\pi^2.
$$

The estimator

$$
\widehat\pi_N=\frac1N\sum_{i=1}^N Y_i
$$

is unbiased and has standard error

$$
\operatorname{SE}(\widehat\pi_N)
=\frac{\sqrt{4\pi-\pi^2}}{\sqrt N}.
$$

### Law of large numbers and central limit theorem

The law of large numbers gives consistency:

$$
\widehat\pi_N\xrightarrow{\mathrm{a.s.}}\pi.
$$

The central limit theorem describes the finite-sample fluctuations:

$$
\frac{\sqrt N(\widehat\pi_N-\pi)}{\sqrt{4\pi-\pi^2}}
\xrightarrow{d}N(0,1).
$$

The first statement says that the estimate approaches the target. The second gives the scale and approximate shape of its remaining error.

### Why the error rate is $N^{-1/2}$

For independent observations with variance $\sigma_Y^2$,

$$
\operatorname{Var}(\overline Y_N)=\frac{\sigma_Y^2}{N}.
$$

Taking a square root gives an error scale proportional to $N^{-1/2}$. On log-log axes,

$$
\log(\mathrm{RMSE})=C-\frac12\log N,
$$

so the theoretical slope is $-1/2$. The notebook estimates RMSE from repeated simulations because the absolute error from a single run need not decrease monotonically.

### Monte Carlo integration

If $U\sim\operatorname{Uniform}(a,b)$, then

$$
E[f(U)]=\frac1{b-a}\int_a^b f(x)\,dx.
$$

Rearranging gives

$$
\int_a^b f(x)\,dx=(b-a)E[f(U)].
$$

The corresponding estimator is

$$
\widehat I_N=\frac{b-a}{N}\sum_{i=1}^N f(U_i).
$$

### Code connection

- `estimate_pi` and `replicate_pi_estimates` implement the geometric estimator.
- `running_mean` displays convergence along one sample stream.
- `root_mean_squared_error` summarizes repeated estimates.
- `estimate_integral` implements uniform-sampling integration and reports an estimated standard error.

## Questions to answer before a v1.0 release

You should be able to answer these without reading the notebook text:

1. Why does independence allow $\operatorname{Var}(S_n)=n$?
2. How can a hitting time be finite almost surely but have infinite expectation?
3. Why does $P(W>t)=P(N(t)=0)$ connect Poisson counts to exponential waits?
4. Why is near-zero empirical correlation weaker than independence?
5. What probability paths are summed by the $(i,j)$ entry of $P^2$?
6. What does $\pi P=\pi$ mean, and why can states still change under stationarity?
7. How does each of $\theta$, $\mu$, and $\sigma$ affect an OU process?
8. Why does Euler--Maruyama use $\sqrt{\Delta t}Z$ for a Brownian increment?
9. What is the difference between the law of large numbers and the central limit theorem?
10. Why does reducing Monte Carlo error by a factor of ten require about one hundred times as many samples?
11. Why are repeated estimates used to measure RMSE?
12. Which conclusions in this repository are theorems, and which are finite-sample observations?

## Reproduction checklist

Before treating the repository as a finished release:

- Run every notebook from a fresh kernel.
- Explain each displayed formula in your own words.
- Change one parameter in every experiment and predict the result before running it.
- Read the corresponding module in `src/` and trace how its outputs enter the notebook.
- Run the full test suite and understand what each test protects.
- Confirm that README values still match the executed notebook outputs.
