# Exact projected-angle kernel and inversion for isotropic Keplerian binaries

## Scope and notation

Let a binary have eccentricity `e in [0,1)`. Its mean anomaly is uniform on a complete orbital period, and its three-dimensional orientation is Haar-uniform on `SO(3)`, independently of phase and eccentricity. Let `Gamma in [0,pi]` be the angle between the two sky-plane projections of its instantaneous separation and relative velocity vectors. Orientations with a zero projected vector have probability zero. Population eccentricities may follow an arbitrary Borel probability measure `F` on `[0,1)`; neither a density nor a parametric form is required.

All statements below concern this ideal phase and orientation sampling. Conditioning on projected separation, proper-motion significance, or any other selection variable need not preserve these assumptions.

Write

\[
j=\sqrt{1-e^2},\qquad q=e/j,\qquad c=|\cos\gamma|.
\]

## 1. A uniform-mixture representation

### Lemma 1 (phase reduction)

In eccentric-anomaly coordinates, up to positive factors irrelevant to direction,

\[
r=(\cos E-e,j\sin E,0),\qquad
v=(-\sin E,j\cos E,0).
\]

Consequently

\[
\frac{r\cdot v}{|r\times v|}=\frac{e\sin E}{j}=q\sin E.
\]

Uniform mean anomaly gives the eccentric-anomaly density

\[
\frac{1-e\cos E}{2\pi}\,dE.
\]

For every bounded measurable function `g` of `sin E`,

\[
\int_0^{2\pi}g(\sin E)\cos E\,dE=0.
\]

This follows by pairing the two preimages of each sine value, or by integrating over the four monotone arcs of sine. Thus any law depending on phase only through `q sin E` can be computed with `E` uniform on `[0,2pi)`.

### Lemma 2 (conditional uniform distribution)

For each phase choose the orthonormal frame whose first vector points along `r` and whose second lies along the component of `v` perpendicular to `r`. Haar-uniform orbital orientation makes this frame Haar-uniform independently of phase. Write its first two vectors as

\[
\widehat r=(\sin\alpha,0,\cos\alpha),\qquad
t=(\cos\alpha\cos\beta,\sin\beta,-\sin\alpha\cos\beta).
\]

The overall sky azimuth has been removed, since it does not affect the measured angle. Here `cos alpha` is uniform on `[-1,1]`, `beta` is uniform on `[0,2pi)`, and the two variables are independent of phase. The velocity direction is proportional to `q sin E rhat + t`. Therefore, exactly,

\[
Y:=\cot\Gamma
=\frac{q\sin E\sin\alpha+\cos\alpha\cos\beta}{|\sin\beta|}.
\]

After Lemma 1, the vector

\[
Z=(\sin\alpha\cos E,\sin\alpha\sin E,\cos\alpha)
\]

is uniform on the unit sphere. At fixed `beta`, the numerator is a linear functional of `Z`, with coefficient-vector norm `sqrt(q^2+cos^2 beta)`. A coordinate of a uniform unit-sphere vector is uniform on `[-1,1]`. Hence

\[
Y\mid\beta,e\sim\operatorname{Unif}[-A_e(\beta),A_e(\beta)],
\qquad
A_e(\beta)=\sqrt{\frac1{j^2\sin^2\beta}-1}.
\tag{1}
\]

The exceptional angles where a denominator or half-width vanishes have zero probability and can be assigned arbitrarily.

Equivalently, if `T=j |sin beta|`, then

\[
f_T(t)=\frac2{\pi\sqrt{j^2-t^2}},\quad 0<t<j,
\qquad
Y\mid T=t\sim\operatorname{Unif}\left[-\frac{\sqrt{1-t^2}}t,
\frac{\sqrt{1-t^2}}t\right].
\tag{2}
\]

This representation also proves normalisation without a separate singular integral calculation.

## 2. The exact likelihood

For `0<e<1`, (2) gives

\[
p_Y(y\mid e)=\frac1\pi
\int_0^{\min(j,(1+y^2)^{-1/2})}
\frac{t\,dt}{\sqrt{1-t^2}\sqrt{j^2-t^2}}.
\]

With `z=sqrt(j^2-t^2)`, the integral has the elementary antiderivative `asinh(z/e)`. Since `asinh(j/e)=acosh(1/e)` and, when `c>e`, `asinh(sqrt(c^2-e^2)/e)=acosh(c/e)`, changing from `Y=cot Gamma` gives

\[
\boxed{
p(\gamma\mid e)=
\frac{\operatorname{acosh}(1/e)
-\mathbf1_{\{|\cos\gamma|>e\}}\operatorname{acosh}(|\cos\gamma|/e)}
{\pi\sin^2\gamma}}
\quad (0<\gamma<\pi).
\tag{3}
\]

The numerically and analytically useful equivalent form is

\[
\boxed{
p(\gamma\mid e)=\frac1{\pi\sin^2\gamma}
\int_{\max(e,|\cos\gamma|)}^1\frac{dt}{\sqrt{t^2-e^2}}.}
\tag{4}
\]

The equality case `c=e` is unambiguous because `acosh(1)=0`. Formula (4) avoids an unnecessary subtraction of two large inverse hyperbolic cosines.

### Circular and endpoint limits

For `e=0`, (4) yields

\[
\boxed{p(\gamma\mid0)=\frac{-\log|\cos\gamma|}{\pi\sin^2\gamma}.}
\tag{5}
\]

Its logarithmic divergence at `gamma=pi/2` is integrable. For each fixed `e<1`,

\[
\lim_{\gamma\downarrow0}p(\gamma\mid e)
=\lim_{\gamma\uparrow\pi}p(\gamma\mid e)
=\frac1{2\pi\sqrt{1-e^2}}.
\]

As `e` tends to one, the law converges weakly to

\[
\tfrac12\delta_0+\tfrac12\delta_\pi.
\]

For example, (1) makes `|Y|` diverge in probability and its sign remains symmetric. An atom at `e=1` can therefore be incorporated as a mathematical closure of the model, although `e=1` itself is not a bound elliptic orbit.

## 3. Arbitrary population mixtures and exact Abel inversion

Let `p_F` be the mixture density for a Borel probability `F` on `[0,1)`. Define, for `0<c<1`,

\[
Y_F(c)=\pi(1-c^2)p_F(\arccos c).
\]

Tonelli's theorem applied to (4) gives

\[
Y_F(c)=\int_c^1 A_F(t)\,dt,
\qquad
A_F(t)=\int_{[0,t)}\frac{F(de)}{\sqrt{t^2-e^2}}.
\tag{6}
\]

The displayed inner integral is finite almost everywhere in `t`; any exceptional null set is immaterial to its locally integrable representative. The strict interval `[0,t)` excludes an atom at `e=t` from the chosen pointwise definition, while an atom can produce an unbounded right-hand limit. For every positive `c`, the integrated quantity in (6) is finite. Hence `Y_F` is locally absolutely continuous on `(0,1)` and

\[
-Y_F'(t)=A_F(t)\quad\text{for almost every }t\in(0,1).
\tag{7}
\]

### Theorem 1 (measure-level inversion and identifiability)

For every `0<r<1`,

\[
\boxed{F([0,r))=
\frac2\pi\int_0^r\frac{t[-Y_F'(t)]}{\sqrt{r^2-t^2}}\,dt.}
\tag{8}
\]

To prove this, substitute (6), apply Tonelli, and use, for every `0<=e<r`,

\[
\int_e^r\frac{t\,dt}{\sqrt{r^2-t^2}\sqrt{t^2-e^2}}=\frac\pi2.
\]

The last identity follows from `u=(t^2-e^2)/(r^2-e^2)`. Thus the entire eccentricity measure is determined by the projected-angle law. The strict endpoint in (8) is deliberate: an atom at `r` is recovered from the right-hand limit of (8) at `r`. The atom at zero is recovered by taking `r` down to zero. No smoothness, parametric-family, or absolute-continuity assumption is used.

Equivalently, if `G` is the distribution of `u=e^2`, then `A_F(sqrt x)` is the classical left Abel integral of `G`:

\[
A_F(\sqrt x)=\int_{[0,x)}\frac{G(du)}{\sqrt{x-u}},
\quad
G([0,x))=\frac1\pi\int_0^x\frac{A_F(\sqrt t)}{\sqrt{x-t}}\,dt.
\]

The theorem extends to a law on `[0,1]` by first recovering its mass at one from the observed atoms at `Gamma=0,pi`, then applying (8) to the remaining subprobability measure.

## 4. Closed power-law likelihood and the unique thermal law

For the usual power-law eccentricity density

\[
f_\alpha(e)=(\alpha+1)e^\alpha,\qquad\alpha>-1,
\]

substitution `e=tu` in (6) gives

\[
A_F(t)=\frac{\alpha+1}{2}
B\!\left(\frac{\alpha+1}2,\frac12\right)t^\alpha.
\]

Consequently

\[
\boxed{p_\alpha(\gamma)=
\frac{B((\alpha+1)/2,1/2)}{2\pi}
\frac{1-|\cos\gamma|^{\alpha+1}}{1-\cos^2\gamma}.}
\tag{9}
\]

In particular,

\[
p_0(\gamma)=\frac1{2(1+|\cos\gamma|)},
\qquad
p_1(\gamma)=\frac1\pi.
\]

The already-known thermal implication `f(e)=2e => uniform Gamma` follows immediately. The new identifiability theorem proves its converse among *all* Borel eccentricity laws under the stated phase and orientation assumptions: a uniform projected-angle law forces `F([0,r))=r^2`, hence `f(e)=2e`.

The limit `alpha` down to minus one in (9) recovers the circular likelihood (5).

## 5. Bounded unbiased observables for every odd moment of circularity

Let

\[
B_k=B(k,1/2),\qquad
W_k(\gamma)=\frac\pi{B_k}
\sin^{2k-1}\gamma\,[1-(2k-1)\cos^2\gamma],\quad k\ge1.
\tag{10}
\]

### Theorem 2 (exact moment identities)

For every `e in [0,1)` and integer `k>=1`,

\[
\boxed{\mathbb E[W_k(\Gamma)\mid e]
=(1-e^2)^{k-1/2}=j^{2k-1}.}
\tag{11}
\]

Proof: put `D_k(c)=c(1-c^2)^(k-1)`. Direct differentiation gives

\[
D_k'(c)=(1-c^2)^{k-2}[1-(2k-1)c^2].
\]

For `k=1` the product on the right simplifies identically to one. Using symmetry of the angle law and (4),

\[
\mathbb E W_k
=\frac2{B_k}\int_0^1Y_e(c)D_k'(c)\,dc
=\frac2{B_k}\int_e^1
\frac{c(1-c^2)^{k-1}}{\sqrt{c^2-e^2}}\,dc.
\]

The boundary terms vanish, including in the circular case because `c log c` tends to zero. Substitution `u=(c^2-e^2)/(1-e^2)` evaluates the last integral as

\[
\frac2{B_k}\cdot\frac{j^{2k-1}}2 B(1/2,k)=j^{2k-1}.
\]

The first observable is particularly simple:

\[
\mathbb E\sqrt{1-e^2}=\frac\pi2\,\mathbb E\sin^3\Gamma.
\tag{12}
\]

Equivalent raw odd-sine moment identities are

\[
\boxed{\pi\,\mathbb E[\sin^{2k+1}\Gamma\mid e]
=B(k,1/2)\sum_{r=0}^{k-1}\frac{j^{2r+1}}{2r+1}.}
\tag{13}
\]

Equation (13) follows from (11), the beta recurrence `B_k/B_(k-1)=2(k-1)/(2k-1)`, and the base case (12).

### Exact boundedness and sampling uncertainty

The maximum of `W_k` is `pi/B_k` at `gamma=pi/2`. For `k>=2`, its minimum is

\[
-\frac{2\pi}{B_k}
\left(\frac{2k-2}{2k+1}\right)^{k+1/2};
\]

for `k=1` its minimum is zero. Since the negative multiplier is less than one, in all cases

\[
\|W_k\|_\infty=\frac\pi{B(k,1/2)}\sim\sqrt{\pi k}.
\tag{14}
\]

To verify the minimum, differentiate `(1-z)^(k-1/2)[1-(2k-1)z]` for `z=cos^2 gamma`; the interior critical point is `z=3/(2k+1)`. Its negative magnitude is below one because, with `r=1-3/(2k+1)`, `2r^((2k+1)/2) < 2 exp(-3/2)<1`.

Thus for independent identically sampled, exact angles, the sample average of `W_k` is an unbiased estimate of the population moment in (11), with

\[
\operatorname{Var}\!\left(\frac1N\sum_iW_k(\Gamma_i)\right)
\le\frac{\pi^2}{N B(k,1/2)^2}\sim\frac{\pi k}{N}.
\]

The observed sample variance gives a sharper estimated standard error. A distribution-free Hoeffding bound follows using the exact range above, or conservatively the range `2pi/B_k`. These statements do not incorporate angle measurement error or selection, both of which must be modelled separately.

### A second proof of identifiability

Let `J=sqrt(1-e^2)` and let `H` be its law. Equation (11) determines every moment of the finite measure

\[
\eta(du)=\sqrt u\,H_{J^2}(du),\qquad0\le u\le1.
\]

Indeed `int u^(k-1) eta(du)=E J^(2k-1)`. The Hausdorff moment problem on `[0,1]` is determinate (polynomials are dense in continuous functions). The measure `eta` determines `H_(J^2)` away from zero, and the total mass determines any remaining atom at zero. This proves population identifiability directly from bounded observables without differentiating an estimated density. Determinacy does not imply that reconstructing a general distribution from finitely many noisy moments is well conditioned.

## 6. Local ill-posedness: what exact inversion does and does not establish

For an eccentricity density written in squared-eccentricity coordinates, `g(u)=f(sqrt u)/(2sqrt u)`, equation (6) says

\[
A_F(\sqrt x)=\int_0^xg(u)(x-u)^{-1/2}\,du.
\]

This is a Riemann--Liouville fractional integral of order one half, up to its normalising constant. Since `A_F(c)=-Y_F'(c)`, inversion entails one ordinary derivative followed by the inverse of this half integral. On compact subintervals of `(0,1)`, where changes of coordinates and the factor `1-c^2` are smooth and nonzero, the overall forward operator therefore has the local smoothing order three halves. Endpoint weights and atoms require separate functional-analytic treatment; the local statement is not a proved global singular-value asymptotic.

Finite-angle samples do not permit direct pointwise differentiation without regularisation. Equations (10)--(14) provide stable estimates of individual moments, with explicit bounds, while the recovery of a fine-grained eccentricity density remains an inverse problem. No unrestricted-noise stability, minimax rate, or end-to-end survey consistency is claimed here.

## 7. A rigorous obstruction to unrestricted inverse stability

### Proposition 3 (failure of continuity in total variation)

Fix `0<delta<1`. For positive integers `n`, define two eccentricity densities

\[
f_n^{\pm}(e)=2e\,[1\pm\delta\cos(2\pi n e^2)],\qquad0<e<1.
\tag{15}
\]

These are nonnegative and integrate to one. Their separation in total variation is independent of `n`:

\[
\|F_n^+-F_n^-\|_{\mathrm{TV}}
=\frac12\int_0^1 4\delta e\,|\cos(2\pi n e^2)|\,de
=\frac{2\delta}{\pi}.
\tag{16}
\]

Nevertheless their projected-angle laws approach one another in total variation. In fact each approaches the uniform projected-angle law generated by the thermal density `2e`.

To see this, write `K(gamma,e)=p(gamma|e)` and set `u=e^2`. The forward densities satisfy

\[
p_n^\pm(\gamma)-\frac1\pi
=\pm\delta\int_0^1K(\gamma,\sqrt u)\cos(2\pi n u)\,du.
\tag{17}
\]

For every interior angle, the kernel as a function of `u` is integrable, since

\[
\int_0^1K(\gamma,\sqrt u)\,du=\frac1\pi
\]

is precisely the thermal identity. The Riemann--Lebesgue lemma makes the integral in (17) tend to zero. Its absolute value is bounded by `1/pi`. Dominated convergence over `gamma in (0,pi)` therefore gives

\[
\|P_n^\pm-\operatorname{Unif}(0,\pi)\|_{\mathrm{TV}}\longrightarrow0,
\qquad
\|P_n^+-P_n^-\|_{\mathrm{TV}}\longrightarrow0.
\tag{18}
\]

In contrast, `||F_n^\pm-F_thermal||_TV=delta/pi`. Thus the inverse is not even continuous in total variation at the uniform angle law when the domain includes all eccentricity densities. Exact population identifiability alone cannot provide unrestricted reconstruction precision. Smoothness, shape constraints, finite-dimensional modelling, or another regularising assumption is needed for uniformly reliable density recovery. This proposition is independent of measurement error and does not rely on a proposed singular-value asymptotic.

## 8. Two precise robustness corollaries

### Known independent angle error

Identify the endpoints of `[0,pi]` to form the circle `R/(pi Z)`; the ideal density is already pi-periodic and even. Suppose an observed angle is

\[
\Theta=\Gamma+\varepsilon\pmod\pi,
\]

where the error law `Q` is known and independent of the true angle and eccentricity. If all Fourier multipliers

\[
\widehat Q(k)=\mathbb E\exp(2ik\varepsilon),\qquad k\in\mathbb Z,
\]

are nonzero, then the eccentricity distribution is still identifiable from the observed-angle population law. Indeed the observed Fourier coefficients are the products `P_hat(k) Q_hat(k)`, so equality of observed laws forces equality of every ideal-angle Fourier coefficient. Fourier uniqueness on the compact circle and Theorem 1 finish the proof.

For example, a Gaussian error of known variance `sigma^2`, wrapped modulo pi, has multiplier `exp(-2k^2 sigma^2)`, which is nonzero for every finite `k`. This establishes exact identifiability; the exponentially small multipliers make recovery of high-frequency features much less stable. No necessity claim is made for the nonzero-multiplier condition. Nor is a physical astrometric-error model asserted to be an independent additive angular error without further justification.

### Known positive selection depending only on the angle

Suppose the probability of retaining a system is a known measurable function `s(gamma)` satisfying `0<s(gamma)<=1` on the support of the ideal angle law. The selected law is

\[
P_s(d\gamma)=\frac{s(\gamma)P(d\gamma)}{\int s\,dP}.
\]

It determines the original law exactly through

\[
P(d\gamma)=
\frac{s(\gamma)^{-1}P_s(d\gamma)}{\int s^{-1}\,dP_s}.
\tag{19}
\]

Therefore known positive angle-only selection preserves eccentricity identifiability. Selection probabilities need not be bounded away from zero for this algebraic conclusion, although very small probabilities can make finite-sample reweighting unstable. Selection on projected separation, speed, inclination, distance, or another correlated variable is outside this corollary unless it has first been reduced to a justified known angle-only probability under the model.

## 9. Necessary observable model checks

### Proposition 4 (moment inequalities)

Write `m_k=E W_k(Gamma)` for `k>=1`. Every projected-angle population satisfying the assumptions obeys

\[
1\ge m_1\ge m_2\ge\cdots\ge0,
\qquad
m_{k+1}^2\le m_km_{k+2},
\]

and, more generally, for every nonnegative integer `r`,

\[
(-1)^r\Delta^r m_k
=\mathbb E[J^{2k-1}(1-J^2)^r]\ge0,
\tag{20}
\]

where `Delta m_k=m_(k+1)-m_k`. Monotonicity and (20) follow directly from `0<=J<=1`; the quadratic inequality follows by Cauchy--Schwarz. These are necessary conditions, not a claim that a finite list characterises the whole model.

The simplest consequence is the sharp observable bound

\[
\boxed{\mathbb E\sin^3\Gamma\le\frac2\pi,}
\tag{21}
\]

with equality if and only if all orbits are circular. Thus excessive concentration of observed angles near `pi/2` cannot be explained solely by changing the eccentricity distribution within the stated model. Measurement errors and selection must be accounted for before applying this diagnostic to data.

For example, for `N` independent exact angles, Hoeffding's inequality gives a distribution-free one-sided specification test: under the model, the probability that

\[
\frac1N\sum_i\sin^3\Gamma_i>
\frac2\pi+\sqrt{\frac{\log(1/a)}{2N}}
\]

is at most `a`, for any `0<a<1`. Other finite differences in (20) give additional tests, with uncertainty derived from their bounded observable combinations.

## Appendix: an exact individual posterior under a flat eccentricity prior

This appendix concerns an explicit prior choice, `f_prior(e)=1` on `[0,1]`; it does not turn one observed projected angle into a prior-free measurement of an individual orbit. Let `c=|cos gamma|<1`. Equations (3) and (9) give

\[
\pi(e\mid\gamma)=\frac2{\pi(1-c)}
\left[\operatorname{acosh}(1/e)
-\mathbf1_{e<c}\operatorname{acosh}(c/e)\right].
\tag{22}
\]

For every real `m>-1`, the posterior moments are

\[
\boxed{\mathbb E[e^m\mid\gamma]
=\frac{B((m+1)/2,1/2)}{\pi(m+1)}
\frac{1-c^{m+1}}{1-c}.}
\tag{23}
\]

Indeed the unnormalised numerator is `p_m(gamma)/(m+1)`, and the denominator is `p_0(gamma)`. In particular,

\[
\mathbb E[e\mid\gamma]=\frac{1+c}{\pi},
\qquad
\mathbb E[e^2\mid\gamma]=\frac{1+c+c^2}{6}.
\]

For a completely elementary posterior CDF, define

\[
L_a(r)=
\begin{cases}
r\operatorname{acosh}(a/r)+a\arcsin(r/a),&0<r<a,\\
a\pi/2,&r\ge a,
\end{cases}
\]

with `L_a(0)=0` and `L_0(r)=0`. The continuous extension at `r=a` is intended. Then, for `0<=r<=1`,

\[
\boxed{\Pr(e\le r\mid\gamma)
=\frac2{\pi(1-c)}[L_1(r)-L_c(r)].}
\tag{24}
\]

This follows from differentiating `r acosh(a/r)+a asin(r/a)`. As `c` tends to one, the CDF tends to `2 asin(r)/pi`.

The posterior mode is `e=c`: the likelihood increases for `e<c` and decreases for `e>c`. Hwang, Ting and Zakamska (2022) already identify this mode relation numerically in their main text. Its exact verification from (3) should be credited as a check of that known observation, not claimed as a new statement. The moment and CDF formulae above remain explicitly conditional on the chosen flat prior.

