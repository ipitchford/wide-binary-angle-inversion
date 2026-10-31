# Novelty-gate report: wide-binary angle identifiability

**Search date:** 28 October 2026.  
**Prepared by:** the literature/novelty agent, with separate mathematical work by the root and projection agents.  
**Corpus boundary:** public web-indexed primary papers and accessible full texts through 28 October 2026, including Oxford Academic, arXiv, DML-CZ, SIAM primary abstracts, publisher abstracts, the complete ADS scan of Tokovinin (1998), and the originating authors' public code repository. Searches used two web search engines and an additional Exa semantic search. This is a bounded search, not a comprehensive bibliographic database or author-confirmed priority investigation.  
**Intended contribution:** an elementary exact likelihood for the projected position–velocity angle of a randomly oriented, equal-time-sampled Kepler orbit; its explicit reduction to a classical Abel transform; and a constructive arbitrary-law injectivity theorem answering the question in Appendix B of Hwang, Ting and Zakamska (2022). Formula-level searches were repeated after the spectral, finite-dual, and final direct-kernel simplifications. The direct likelihood and Abel bridge are now the main contribution; the older spectral construction is retained below only as an audit trail.

**Gate result:** no exact likelihood, Abel bridge, or complete arbitrary-law injectivity collision was found in the inspected corpus. The result is **bounded uncertainty, with a strong basis to proceed**. The original observational idea, the thermal example, numerical inversion, standard spherical geometry, and Abel inversion are established prior art. No categorical claim of being the first proof is justified by this search alone.

## 1. Frozen claim and exclusions

Let an arbitrary Borel probability law on eccentricity $e\in[0,1)$ generate a population of bound Kepler orbits. Conditional on eccentricity, mean anomaly is uniform on $[0,2\pi)$, and orbital orientation is independent Haar measure on $SO(3)$. The observable $\Gamma\in[0,\pi]$ is the angle between the two sky-plane projections of relative position and relative velocity. The target is that the law of $\Gamma$ uniquely determines the eccentricity law, including atoms.

The intended statement is distribution-level. It does not identify the eccentricity of an individual binary from one angle, establish good finite-sample resolution at every eccentricity, or validate the assumptions for a selected Gaia catalogue. Numerical rank of a discretised kernel and identifiability of a single power-law family would not resolve the target.

Permitted priority language is “we give a constructive solution to the identifiability question posed in Hwang et al. (2022)”, once the mathematical proof is verified, together with a precise description of the inspected prior work. The word “first” is not warranted by negative search evidence. A formula's algebraic simplicity does not by itself make its mathematical content new.

## 2. Originating question and inspected astronomy sources

### Hwang, Ting and Zakamska (2022)

**Primary reference:** H.-C. Hwang, Y.-S. Ting and N. L. Zakamska, *The eccentricity distribution of wide binaries and their individual measurements*, MNRAS **512**, 3383–3399, DOI [10.1093/mnras/stac675](https://doi.org/10.1093/mnras/stac675). Full author manuscript: [arXiv:2111.01789v2](https://arxiv.org/pdf/2111.01789).

**Inspected passages:** Sections 2.1–2.5; Appendix A; Appendix B (PDF pp. 18–19). Section 2 samples mean anomaly uniformly and takes random orbital orientations. Appendix A writes the projected likelihood as a three-dimensional integral over true anomaly and orientation. Appendix B explicitly asks whether two distinct eccentricity distributions can give the same observed angle distribution. It demonstrates distinction numerically for the power-law family and leaves the general question unresolved there. The same appendix identifies the loss of uniqueness when orientation is not random.

This is direct evidence that the question was unresolved in that paper. It is not proof that nobody answered it later. The 2022 paper also credits Tokovinin (1998) for the original observational idea.

**Known ingredient that must receive credit:** Eqs. (11)–(16) already specify the intrinsic three-dimensional angle density by summing its two true-anomaly branches. The compact kernel obtained in the present work,

$$p(\mu\mid e)=\frac{\sqrt{1-e^2}}{\pi(1-\mu^2)\sqrt{e^2-\mu^2}},\qquad |\mu|<e,$$

with $\mu=\cos\theta$, is a simplification of that same known intrinsic forward density. Its existence or its physical interpretation must not be presented as a newly discovered forward law. In contrast, Section 2.5 explicitly treats an analytic expression for the **projected** likelihood as challenging and evaluates it by simulation; Appendix A gives a multidimensional integral evaluated numerically. The proposed contribution is the closed projected likelihood, its reduction to an Abel inverse, and the resulting arbitrary-law injectivity and observable moments.

The authors' [public code repository](https://github.com/HC-Hwang/Eccentricity-of-wide-binaries), inspected through its README, supplies precomputed $p(\gamma\mid e)$ grids generated by `binary_sampling`, at eccentricities from 0 to 0.99 in steps of 0.01. This supplies a concrete computational comparison: the exact formula evaluates the idealised likelihood directly without that simulation grid. It does not by itself establish a measured runtime advantage or remove observational modelling requirements.

### Tokovinin (1998), and Tokovinin–Kiyaeva (2016)

**Primary references:** A. A. Tokovinin, *On the Distribution of Orbital Eccentricities for Wide Visual Binary Stars*, Astronomy Letters **24**, 178–179 (1998), with Russian original in *Pis'ma v Astronomicheskii Zhurnal* **24**, 217–219; [publisher abstract](https://www.pleiades.online/cgi-perl/search.pl?name=letters&number=2&page=178&type=abstract&year=98). A. Tokovinin and O. Kiyaeva, *Eccentricity distribution of wide binaries*, MNRAS **456**, 2070–2079 (2016), DOI [10.1093/mnras/stv2825](https://doi.org/10.1093/mnras/stv2825).

The complete two-page English paper was ultimately obtained from the [primary ADS scan](https://articles.adsabs.harvard.edu/pdf/1998AstL...24..178T). Both pages 178–179 were visually inspected independently by the stereology/proof agent and the literature agent. This closes the main historical access gap in the earlier draft of this report.

Page 178 proposes the angle method and displays Monte Carlo curves from 10,000 realisations at eccentricities 0.2, 0.5, 0.7 and 0.9. Page 179 describes fitting a finite linear combination of those curves to an observed angle histogram. It uses the fitted weights to suggest the thermal law, then verifies the approximately flat angle distribution by Monte Carlo. Its only numbered formula is the thermal law $f(e)de=2e\,de$. It gives no elementary projected likelihood, no Abel bridge, and no general exact injectivity proof. It also cautions that the angle law is insensitive to small eccentricities at the available precision. The observational idea and finite-mixture inference must be credited to this paper; theoretical identifiability does not erase its practical sensitivity concern.

The 2016 paper's method and appendix, inspected through full-text search extracts, reconstruct a binned eccentricity distribution by comparing joint angle–speed histograms with simulations and fitting finite bin weights. It is established prior art for numerical inference. It does not supply a general analytic injectivity theorem in the inspected passages.

### Tokovinin (2020)

**Reference:** *Eccentricity distribution of wide low-mass binaries*, MNRAS **496**, 987–993, DOI [10.1093/mnras/staa1639](https://doi.org/10.1093/mnras/staa1639), [author preprint](https://arxiv.org/abs/2004.06570).

Inspected method passages use the projected angle together with normalised speed. The thermal law gives a flat angle distribution. This established forward implication is not a new claim. A proof that the thermal law is the unique law yielding flat projected angles follows from the proposed general theorem, but that converse should still be searched separately before any independent priority claim.

### Benisty, Evans and Davis (2023)

**Reference:** *The Eccentricity Distribution of Wide Binaries*, MNRAS Letters **518**, L51–L56, DOI [10.1093/mnrasl/slac134](https://doi.org/10.1093/mnrasl/slac134), [arXiv:2210.06429](https://arxiv.org/abs/2210.06429).

Section 3, inspected in full-text search results, derives an Abel-type integral for the distribution of instantaneous **three-dimensional speed normalised by circular speed**, together with acceleration statistics. This is a relevant neighbouring inverse problem, but it uses a different observable from the projected angle alone. The ordinary Abel machinery is classical and should not be claimed as novel.

### Makarov (2025): later status check

**Reference:** V. V. Makarov, *Distributions of Wide Binary Stars in Theory and in Gaia Data. I. Generalized Ambartsumian (1937) Approach and the Family of Power-law Distributions of Eccentricity*, AJ **169**, 113, [arXiv:2501.02587](https://arxiv.org/html/2501.02587v1).

Sections VI–VII were inspected. This paper describes the motion angle as depending on the intrinsic eccentricity distribution and uses forward Monte Carlo simulation with a prescribed power-law family and measurement-error model. It does not give the general analytic inverse in those sections. Its Part II, [arXiv:2507.19273](https://arxiv.org/html/2507.19273v1), reconstructs semimajor-axis distributions from apparent separations and addresses a distinct inverse problem. These later papers support the continuing relevance of angle-based inference; they do not establish the absence of another later solution.

### Brandt (2026): contemporary closed-form orbital statistics

**Reference:** T. D. Brandt, *Closed-Form Statistical Relations Between Projected Separation, Semimajor Axis, Companion Mass, and Host Acceleration*, [arXiv:2601.14688v2](https://arxiv.org/html/2601.14688v2), revision dated 15 April 2026. The full HTML was inspected.

This paper derives analytic distributions for mass–acceleration–separation relations and an elliptic-integral density for projected separation divided by semimajor axis. These are distinct observables from the projected position–velocity angle. Its Section 5 still integrates the fixed-eccentricity likelihood numerically over a general eccentricity law. No projected angle likelihood or angle-to-eccentricity-law inverse appears in the inspected paper. It provides useful contemporary context for replacing Monte Carlo orbital calculations by analytic densities, not a novelty certificate for the present result.

### Savransky, Cady and Kasdin (2011): standard geometry and neighbouring marginals

**Reference:** D. Savransky, E. Cady and N. J. Kasdin, *Parameter Distributions of Keplerian Orbits*, ApJ **728**, 66, DOI [10.1088/0004-637X/728/1/66](https://doi.org/10.1088/0004-637X/728/1/66), [full manuscript](https://arxiv.org/pdf/1012.1526). Sections 2–4 were inspected.

This paper derives phase, radius and apparent-separation distributions for arbitrary eccentricity and semimajor-axis laws, and establishes isotropy of the instantaneous separation direction. It does not derive the position–velocity projected angle law in those sections. The uniform distribution of a component of a uniform spherical vector, used in the new short proof, is standard geometry and receives no originality claim.

## 3. Definitions and translation checks

| Literature notation or construction | Present normalisation | Check |
|---|---|---|
| Hwang $\gamma_{3D}$, $\Gamma_{3D}=\cos\gamma_{3D}$ | $\theta$, $\mu=\cos\theta$ | Same intrinsic angle and cosine. |
| Hwang $\gamma\in[0,\pi]$ | $\Gamma$ | Same projected unsigned angle. |
| Tokovinin folded $\gamma\in[0,\pi/2]$ | $\min(\Gamma,\pi-\Gamma)$ | Equivalent for the time-reversal-symmetric Kepler population; sine moments are unchanged. |
| Uniform mean anomaly | $dP(E)=(1-e\cos E)dE/(2\pi)$ | Equal-time phase sampling, not uniform true or eccentric anomaly. |
| Isotropic orbital normals with uniform pericentre | Haar $SO(3)$ | Includes both random plane and random orientation within that plane. |
| $a_{2n}=E[P_{2n}(\mu)]$ | Intrinsic even Legendre moments | Odd moments vanish for the Kepler phase mixture. |
| $b_{2k}=E[\cos(2k\Gamma)]$ | Projected even Fourier moments | Specifies the symmetric observed angle law. |

## 4. Formula fingerprints and rerun

### Main fingerprint: elementary projected likelihood

The decisive structural simplification is the following expression, for $0<e<1$ and $0<\gamma<\pi$:

$$p(\gamma\mid e)=\frac{1}{\pi\sin^2\gamma}
\int_{\max(e,|\cos\gamma|)}^1\frac{dt}{\sqrt{t^2-e^2}}
=\frac{\operatorname{acosh}(1/e)-\mathbf1_{\{|\cos\gamma|>e\}}\operatorname{acosh}(|\cos\gamma|/e)}{\pi\sin^2\gamma}.$$

Writing $j=\sqrt{1-e^2}$ and $s=\sin\gamma$ gives the equivalent fingerprint

$$p(\gamma\mid e)=\frac{\operatorname{atanh}j-\mathbf1_{\{s<j\}}\operatorname{atanh}\sqrt{\frac{j^2-s^2}{1-s^2}}}{\pi s^2}.$$

The circular limit is

$$p(\gamma\mid0)=\frac{-\log|\cos\gamma|}{\pi\sin^2\gamma}.$$

The equivalent logarithmic numerator for $c=|\cos\gamma|>e$ is
$\log[(1+\sqrt{1-e^2})/(c+\sqrt{c^2-e^2})]$; for $c\le e$ it is
$\log[(1+\sqrt{1-e^2})/e]$. These are alternative spellings of one object, not distinct discoveries.

The direct proof's probabilistic fingerprint is a uniform scale mixture: conditional on $\beta$ uniform on $[0,2\pi)$,

$$\cot\Gamma\sim \operatorname{Unif}[-A_\beta,A_\beta],\qquad
A_\beta=\sqrt{\frac{1}{(1-e^2)\sin^2\beta}-1}.$$

This representation follows after pairing equal-sine eccentric-anomaly phases so that the unequal-time weight cancels, then recognising a linear functional of a uniform spherical vector. Uniform scale mixtures and spherical coordinate marginals are classical. The candidate new content is their exact connection to this Keplerian projected observable, not a general theorem about those classical probability families.

### Main fingerprint: Abel bridge and power-law family

For an eccentricity density $f$, define $h(c)=\pi(1-c^2)p_\Gamma(\arccos c)$, $0<c<1$. The elementary kernel implies

$$-h'(c)=\int_0^c\frac{f(e)}{\sqrt{c^2-e^2}}\,de.$$

This becomes a classical left Abel integral after the squared-coordinate change. For arbitrary Borel eccentricity laws the same bridge is interpreted almost everywhere or through nonnegative integrals; this report does not substitute formal differentiation for the separate measure-theoretic proof. The original target is uniqueness for arbitrary laws, not only smooth densities. The genuinely substantive bridge is from an orbital observable averaged over phase and viewing angle to an explicitly invertible one-dimensional integral. Abel inversion itself is established machinery.

For the exact family used in Hwang et al., $f(e)=(\alpha+1)e^\alpha$, $\alpha>-1$, the resulting fingerprint is

$$p_\Gamma(\gamma)=\frac{B((\alpha+1)/2,1/2)}{2\pi}
\frac{1-|\cos\gamma|^{\alpha+1}}{\sin^2\gamma}.$$

Important exact checks are

$$\alpha=0:\quad p_\Gamma(\gamma)=\frac{1}{2(1+|\cos\gamma|)},\qquad
\alpha=1:\quad p_\Gamma(\gamma)=\frac1\pi.$$

The thermal second case is known. The general power-law formula and uniform-eccentricity first case were included in the formula search; no matching expression was found in the inspected corpus. A corresponding moment fingerprint is

$$\pi E[\sin^{2k+1}\Gamma\mid e]
=B(k,1/2)\sum_{r=0}^{k-1}\frac{(1-e^2)^{r+1/2}}{2r+1},\qquad k\ge1.$$

The final simplification expresses these as single bounded observable functions:

$$W_k(\gamma)=\frac\pi{B(k,1/2)}\sin^{2k-1}\gamma\,[1-(2k-1)\cos^2\gamma],
\qquad E[W_k(\Gamma)\mid e]=(1-e^2)^{k-1/2}.$$

A last formula/terminology rerun for the sine-cubed mean, circularity moments and unbiased angle observables found no matching Keplerian identity in the returned results. The statistical facts that a sample mean of an unbiased bounded observable is unbiased and admits elementary variance bounds are standard; the candidate contribution is the explicit observable and its orbital expectation.

The acosh, atanh, logarithmic, circular, uniform-eccentricity, power-law and cotangent-mixture forms were searched after this simplification. Searches also used the aliases “projected flight-path angle”, “motion angle”, “relative motion direction”, “projected orthogonal vectors”, and “uniform scale mixture”. The closest retrieved astronomy expressions involving atanh were Benisty et al.'s normalised **speed** distributions. Their observable and transform were checked as different. No exact projected-likelihood or Abel-bridge collision was located.

### Audit history: forward coefficient matrix

The following spectral fingerprints belong to an earlier, longer proof route. They are retained for reproducibility of the gate and are not required by the direct likelihood proof.

The proposed projection transform has

$$\frac{b_{2k}}{k^2}=\sum_{n\ge k}s_n a_{2n}R_{n-k}R_{n+k},\quad
s_n=\frac{(4n+1)\Gamma(n)^2}{4\Gamma(n+3/2)^2},\quad
R_j=\frac{\Gamma(j+1/2)}{\Gamma(j+1)}.$$

Exact columns for $n=1,2,3,4$ are

$$\left(\frac56\right),\qquad
\left(\frac1{10},\frac7{10}\right),\qquad
\left(\frac{13}{420},\frac{26}{175},\frac{429}{700}\right),\qquad
\left(\frac{17}{1260},\frac{187}{3150},\frac{2431}{14700},\frac{2431}{4410}\right).$$

The entries below the diagonal vanish. Every displayed entry is positive. Each finite truncation is upper triangular with nonzero diagonal. Its determinant is the product of those diagonal entries; this is a fingerprint, not the injectivity proof for the infinite transform. Three exact finite truncations are recorded in `target.yaml`.

The first six $\pi s_n$ values are

$$\frac{20}{9},\ \frac{16}{25},\ \frac{3328}{11025},\ \frac{17408}{99225},\ \frac{65536}{571725},\ \frac{6553600}{81162081}.$$

The first row of the matrix is

$$\frac56,\ \frac1{10},\ \frac{13}{420},\ \frac{17}{1260},\ \frac7{990},\ \frac{25}{6006},\ldots,$$

and equals $4(4n+1)/[(2n-1)(2n)(2n+1)(2n+2)]$. The generating function used to recognise the Fourier connection is the classical Legendre generating function $(1-2xt+t^2)^{-1/2}$. No novelty is attached to that identity.

### Audit history: finite dual functions

The simplification that triggered a second gate was the existence of a finite test function $\psi_n$ satisfying $E[\psi_n(\Gamma)\mid\mu]=P_{2n}(\mu)$, with the first three functions

$$\begin{aligned}
\psi_1(\gamma)&=1-\frac{3\pi}{4}\sin^3\gamma,\\
\psi_2(\gamma)&=1-\frac{75\pi}{16}\sin^3\gamma+\frac{315\pi}{64}\sin^5\gamma,\\
\psi_3(\gamma)&=1-\frac{3675\pi}{256}\sin^3\gamma+\frac{19845\pi}{512}\sin^5\gamma-\frac{51975\pi}{2048}\sin^7\gamma.
\end{aligned}$$

Equivalent general fingerprints are

$$\psi_n=1-\frac{8\Gamma(n+3/2)^2}{\Gamma(n)^2}H_n,\qquad
H_n''(\gamma)=\sin\gamma\,P_{2n}(\cos\gamma),\quad H_n(0)=H_n'(0)=0,$$

and the terminating hypergeometric form supplied by the proof agent,

$$H_n(\gamma)=\frac{\sin^3\gamma}{6}\,{}_3F_2\!\left(1-n,n+\frac32,\frac32;2,\frac52;\sin^2\gamma\right).$$

The especially simple eccentricity identity is

$$E\sqrt{1-e^2}=\frac\pi2 E\sin^3\Gamma.$$

The weighted moments reduce to

$$E[(1-\mu^2)\mu^{2n}]=\frac{\binom{2n}{n}}{4^n}E[\sqrt{1-e^2}\,e^{2n}].$$

The central-binomial factor is the classical arcsine moment sequence. Its presence is an ingredient identity, not evidence of novelty. The proposed content lies in connecting observable angle moments to all moments of the weighted eccentricity law.

## 5. Cross-disciplinary searches and verified distinctions

### Dihedral-angle stereology

**Inspected primary article:** V. Horálek, *Stereology of dihedral angles*, Applications of Mathematics **45** (2000), 411–417, DOI [10.1023/A:1022329432636](https://doi.org/10.1023/A:1022329432636), [full text](https://www.dml.cz/bitstream/handle/10338.dmlcz/134449/AplMat_45-2000-6_2.pdf).

**Related primary sources located:** J. A. Reeds and J. P. Butler, *Stereology of Dihedral Angles II: Distribution Function*, SIAM J. Appl. Math. **47** (1987), 678–687, DOI [10.1137/0147046](https://doi.org/10.1137/0147046); J. P. Butler and J. A. Reeds, Part I, 670–677, DOI [10.1137/0147045](https://doi.org/10.1137/0147045); R. E. Miles, *Dihedral angle distributions*, Acta Stereologica **6** (1987), 19–24. The Reeds–Butler primary abstracts were inspected; their full text and Miles's full text remained unavailable through these retrieval paths at the time of this draft.

There is a real geometrical alias: the normals to the traces of two planes in a section plane are the projections of the original plane normals. However, the probability law of the section matters. Horálek pp. 413–414 explicitly conditions on the random section intersecting an edge element and obtains orientation density proportional to $\sin\theta\cos\theta$. The astronomical pair is viewed with unweighted isotropic orientation, proportional to $\sin\theta$. Thus the inspected stereological kernel is not the same probability operator. Simply naming both “random sections/projections” would obscure that difference. No exact collision with the stated unweighted astronomy inverse was established.

Horálek reports that Miles formulates intrinsic cot-angle density as a solution of a symmetric first-kind Fredholm equation. Because that primary full text was not inspected, its complete scope remains a priority-audit limitation. It should be obtained before journal submission.

This literature is now peripheral to the primary direct Kepler likelihood and Abel bridge: the short proof does not assert a new inverse for every three-dimensional pair-angle distribution. The missing Miles full text should still constrain any broader stereological priority claims.

The separate stereology audit also inspected M. A. Fortes, *A general method to solve the integral equations of stereology*, Acta Stereologica **10**(2) (1991), 193–202, [full primary PDF](https://www.ias-iss.org/acta/ActaStereologica_pdf/AS_2008.pdf). Pages 196–198 discretise the unknown distribution into classes and solve a finite system $Y=PX$; p. 198 assumes that $P$ is nonsingular. This is prior art for numerical inversion, not the arbitrary-measure identifiability argument required here. Fortes cites M. A. Fortes and A. C. Ferro (1988), Metallurgical Transactions A **19**, 1147–1151, DOI [10.1007/BF02662574](https://doi.org/10.1007/BF02662574), for a dihedral application; that cited article was not separately inspected in full.

### Fixed-axis orientation inverse

**Inspected primary reference:** P. A. Kloza and J. A. Elliott, *Inverse projection of axisymmetric orientation distributions*, EPL **143**(2) (2023), **22002**, DOI [10.1209/0295-5075/ace4ee](https://doi.org/10.1209/0295-5075/ace4ee), [full publisher text](https://academic.oup.com/epl/article/143/2/ace4ee/8765552). The [original published PDF in the authors' institutional repository](https://api.repository.cam.ac.uk/server/api/core/bitstreams/65846b75-6b97-4840-af9b-4df6a1bce5b7/content) explicitly gives article number 22002. The `ace4ee` string is the DOI suffix and the current publisher's URL label, not a conflicting article number.

This paper gives an Abel inverse and a Chebyshev–Legendre matrix representation. It explicitly assumes that the projection plane contains the symmetry axis. The present model averages the symmetry axis itself over orientation. The difference is verified at the coefficient level: Kloza–Elliott Eq. (23) has coefficients $5/4,-3/8,13/64,\ldots$ for $E[T_2]$ as a function of intrinsic Legendre moments; the present random-pair operator has $5/6,1/10,13/420,\ldots$. The matrices are not the same under the normalisations checked here.

Kloza–Elliott identify the classical Leadbetter–Norris and Deutsch projection/inversion literature. Those sources are appropriate conceptual background, but the present bounded search has not established that their fixed-axis inverse yields the full random-axis identity by a known composition. A claim of novelty for all angular inversion would therefore be wrong.

### Other aliases

The search also covered angular correlations of isotropic vector pairs, marginal azimuth differences, spherical random projection, and projected orientation distributions in fibre/materials analysis. Much “spherical random projection” literature uses Gaussian random matrices, not orthogonal projection onto an isotropically chosen plane in a fixed three-dimensional space. This different probability model must be checked before transferring an inverse or a kernel.

## 6. Reproducible query log

All queries below were run on 28 October 2026. Quotation marks indicate literal phrases where present. Negative results mean no relevant collision appeared in the returned results, not that no such paper exists.

| Family | Exact queries or representative strings actually run | Outcome inspected |
|---|---|---|
| Originating astronomy question | `wide binary eccentricity distribution projected angle inverse kernel identifiability analytic`; `"eccentricity distribution" "APPENDIX B"`; `"eccentricity distribution" "unique mapping"`; `"eccentricity distribution" "injective"` | Hwang Appendix B explicit question; no later general proof returned. |
| Historical astronomy | `1998 "eccentricity distribution" "angle" wide visual binaries Tokovinin`; `"On the Distribution of Orbital Eccentricities" pdf`; `"1998AstL...24..178T"`; `"Tokovinin" "1998" "178" pdf eccentricities` | 1998 publisher abstract and 2016 numerical method; 1998 scan initially inaccessible, ultimately retrieved and both pages visually inspected. |
| Russian historical check | `Токовинин распределении эксцентриситетов широких визуально двойных 1998`; `"Токовинин" "1998" "217"` | Bibliographic confirmation, not a complete 1998 text. |
| Later astronomy status | `wide binary eccentricity distribution inversion uniqueness proof 2026`; `wide binary velocity separation angle 2025 eccentricities`; `wide binary eccentricity analytical angle kernel arxiv`; `Hwang Ting Zakamska eccentricity distribution citing analytic 2024 2025 2026` | Makarov 2025 I–II, Benisty et al. 2023, Galactic-evolution papers; no general angle inverse found. Recent searches used a date window. |
| Nonparametric/Abel aliases | `"eccentricity distribution" "nonparametric" "wide"`; `"eccentricity distribution" "Abel transform"`; `"projected" "angle" "eccentricity" "Abel"` | Benisty speed kernel; no identical projected-angle theorem. |
| Random-angle geometry | `"projected angle" "inverse" "isotropic" distribution`; `"distribution of projected angles"`; `"random projections" "angle distribution" inversion`; `"angle distributions" "stereology"`; `"angle distribution" "injectivity" projection` | Dihedral stereology and different fixed-axis projection problems. |
| Stereological sources | `"Stereology of dihedral angles" "inverse"`; `"Stereology of dihedral angles" "Part II" pdf`; `"Dihedral angle distributions" Miles pdf`; `"Miles" "Dihedral angle distributions"` | Horálek full text, Reeds–Butler abstracts, Miles citation. |
| Spectral simplification | `"projected angle" "Legendre"`; `"angle distribution" "Legendre" "inverse" projection`; `"azimuth difference" "Legendre"`; `"azimuthal correlation" "Legendre" "projection"` | Kloza–Elliott 2023, whose fixed-axis geometry and matrix were verified as different. |
| Exact matrix entries | `"13/420" "26/175"`; `"5/6" "1/10" "7/10" "Legendre"`; `"13/420" "17/1260"`; `"projected" "angle" "5/6" "7/10"` | No matching angular operator returned. |
| Gamma ratios | `"angle" "Gamma(n)" "Gamma(n+3/2)"`; `"Gamma(n)" "Gamma(n+3/2)" "Legendre"`; `"projection" "Gamma" "Legendre" "angle" "inverse"` | Standard special-function/background identities; no exact inverse collision. |
| Finite-dual rerun | `"75/16" "315/64"`; `"sin^3" "projected angle"`; `"sin^3" "dihedral" stereology`; `mean sine cubed projected angle random projection vectors`; `"projected" "angle" "sine cubed"` | No matching finite dual or sine-cubed identity returned. |
| Eccentricity-moment rerun | `"eccentricity" "sin^3" "angle" distribution`; `"eccentricity" "sqrt(1-e^2)" "angle" average`; `"eccentricity" "sin" "cube" "angle"` | Standard orbital averaging and unrelated uses of sine-cubed inclination; no matching projected-angle moment returned. |
| Elementary projected-kernel rerun | `"eccentricity" "angle" "arccosh"`; `"binary" "angle" "eccentricity" "arccosh"`; `"eccentricity" "arctanh" "probability"`; `"eccentricity" "arccosh" "likelihood"`; `"eccentricity" "angle" "arcosh"`; `"projected angle" "atanh"` | No exact likelihood collision; hyperbolic orbital-anomaly relations were excluded as a different object. Benisty speed atanh formula was checked as a different observable. |
| Logarithmic and circular fingerprints | `"projected" "angle" "log" "sin" "circular orbit"`; `"projected angle" "log" "cos" "sin" distribution`; `"orthogonal vectors" "angle" "log" projection distribution`; `random orthogonal vectors projected angle probability density logarithm`; `"eccentricity distribution" "logarithm"` | No matching general or circular projected likelihood returned. The circular case alone should not receive a categorical novelty claim. |
| Power-law and uniform fingerprints | `"eccentricity" "1+|cos"`; `"eccentricity" "1+cos" "angle" distribution`; `"eccentricity" "angle distribution" "cos" "sin" analytic`; `"binary" "angle distribution" "analytic" eccentricity` | Hwang simulations and unrelated formulae; no matching closed power-law law located. |
| Cotangent mixture and alternative names | `"eccentricity" "projected" "cot" distribution angle`; `"projected flight-path angle" distribution`; `"eccentricity" "mixture" "uniform" "cot"`; `"angle" "scale mixtures of uniforms"`; `"angle distribution" "arcsine" "Kepler"`; `"motion angle" "distribution" "analytic"` | No exact Kepler uniform-mixture representation returned. Generic mixture machinery is not claimed as new. |
| Latest rerun | `"eccentricity" "angle" "analytic" "wide binaries"`; `"eccentricity distribution" "inversion" "angle"`; `"eccentricity" "acosh"`, with 730-day window | No resolving paper returned in the bounded date-filtered search; publication/crawl dates were distinguished. |
| Semantic search | “A research paper that analytically derives the probability density of the projected angle between relative velocity and separation in randomly oriented Keplerian binary orbits, using an inverse hyperbolic cosine or logarithm, or an Abel integral inverse of the eccentricity distribution.” | Exa located Brandt 2026, Makarov 2025, Savransky 2011, Benisty 2023 and Hwang 2022. Brandt's complete text was inspected and its distinct observables verified. |
| Public implementation | `site:github.com "eccentricity" "arccosh" "gamma"`; `site:github.com "eccentricity" "arctanh" "angle"`; `site:github.com "Eccentricity-of-wide-binaries"` | Hwang repository README describes precomputed simulation grid; no indexed matching exact kernel implementation returned. |
| Bounded circularity observables | `"eccentricity" "sin^3" "pi/2"`; `"circularity" "projected angle" "moment"`; `"eccentricity" "sin^{2k"`; `"eccentricity" "unbiased" "angle" "moment"` | No matching orbital moment observable returned; unrelated ellipticity and conic-section formulae excluded. |

## 7. Missingness audit

1. The original Tokovinin (1998) access gap is **closed**: both pages of the primary ADS scan were independently visually inspected by two agents. The paper contains a numerical finite-mixture approach and no exact projected kernel or arbitrary-law injectivity theorem. Earlier source-note statements to the contrary refer to retrieval difficulty at that earlier time, not to the final audit status.
2. Miles (1987) and Reeds–Butler (1987) full texts were not fully inspected; the latter's primary abstracts were read and Horálek supplies explicit sampling conventions and references. This remains a limitation on broad random-angle inversion priority claims. The final direct Kepler composite theorem does not require such a broad claim.
3. Search engines do not index mathematical notation reliably. Literal fraction and gamma-function queries have low recall. Formula-level negative evidence is therefore weaker than a verified formula mismatch.
4. No paid MathSciNet, Web of Science, Scopus or exhaustive ADS citation export was available. Citation searches were web based. No systematic search of all theses, proceedings or all non-English literature was completed.
5. No author enquiry was sent. Author confirmation and expert bibliographic review would remain valuable for priority, particularly for historical angular-projection mathematics.
6. The current search is bounded by public indexing as retrieved on the search date. It cannot certify the nonexistence of an unindexed 2026 manuscript.

## 8. Outcome by contribution unit

| Unit | Outcome | Evidence and allowed language |
|---|---|---|
| Angle-based eccentricity inference as an observational idea | **COLLISION** | Proposed by Tokovinin (1998), developed numerically in 2016/2020/2022. Credit those papers. |
| Intrinsic three-dimensional angle forward density | **COLLISION / algebraic simplification** | Hwang Eqs. (11)–(16) already define the same density. The compact expression may be used as a lemma with attribution. |
| Thermal law gives flat projected angles; modal relation $|\cos\gamma|=e$ | **COLLISION as forward statements** | The thermal implication is historical; Hwang numerically identifies the modal relation. An exact proof may be a useful corollary but not a newly observed statement. |
| Classical Abel and Hausdorff moment ingredients | **COLLISION / standard machinery** | No novelty claim for the underlying transforms or determinacy theorem. |
| Elementary projected likelihood and direct Kepler-to-Abel bridge | **BOUNDED UNCERTAINTY; strong basis to proceed** | No exact collision in the inspected astronomy sources or formula reruns. Hwang evaluates the same projected likelihood by simulations and multidimensional integration; Tokovinin 1998 full text is numerical. Describe the exact formula and bridge, with no categorical priority wording. |
| Full arbitrary-law projected-angle injectivity | **BOUNDED UNCERTAINTY; clear to investigate within inspected corpus** | A complete proof was not found in the bounded search. The source explicitly posed the question in 2022. A rigorous new proof can be described as answering that paper's question, with historical priority reserved. |
| Explicit projected law for all eccentricity power laws | **BOUNDED UNCERTAINTY** | Same parametric family is known and widely used; no closed angle-law expression was located. Credit the family and thermal case; distinguish analytic evaluation from inventing the model. |
| Finite sine-polynomial dual functions and weighted eccentricity moments | **BOUNDED UNCERTAINTY** | Exact fingerprints were searched; no collision was found. Related fixed-axis and section-biased kernels were distinguished. Avoid “first-ever” language. |
| Reproducibility and exact verification artefacts | **CLEAR** | Independently generated work in this session, subject to the mathematical checks described elsewhere. |

## 9. Decision

**Proceed with the exact-likelihood paper and independent proof verification, while retaining bounded priority language.** The originating question is documented directly in the 2022 paper. The complete 1998 source supplies the original idea and numerical inverse, not the claimed exact formula or theorem. Contemporary papers inspected through 2026 address neighbouring observables or continue angle-based forward simulation. The strongest contribution unit is the unexpectedly elementary Kepler likelihood together with the explicit Abel bridge and arbitrary-law uniqueness proof; merely renaming Abel inversion would not be a contribution.

This bounded gate supports the claim that the manuscript gives a mathematically substantive answer to the question posed by Hwang et al. It does not certify historical priority. Before external submission, a domain expert should check novelty, practical value, and the treatment of selection and noise. The residual old stereology gaps matter principally if the manuscript separately claims a new universal projection inverse.

## 10. Rerun triggers

Rerun the gate if a further transformation identifies the **composite Kepler likelihood** with a known named distribution, if a new correction or theorem materially changes the claim, or if a later paper directly cites and resolves Hwang Appendix B. The gate was rerun for the gamma-ratio spectral matrix, the finite sine-polynomial duals, and finally the elementary acosh/atanh likelihood, uniform-scale-mixture representation, Abel bridge and explicit power-law family. Identifying the inverse step as classical Abel inversion has already been incorporated into the contribution accounting.


## Addendum, 31 October 2026: results added in revision

The revision adds real-order moment observables, a dual-observable uniqueness
theorem, Fourier formulas for even eccentricity moments, power-law index
estimators with Fisher information, a weak-continuity result, a thermal
selection-invariance proposition and a separation-binning proposition.

**Checks performed.** Five internal reviewers searched the astronomy and
inverse-problem literature on 31 October 2026 (web, arXiv, Crossref, Semantic
Scholar). Added comparisons: Benisty, Evans and Davis (2022), analytic speed and
acceleration laws for the same eccentricity problem (different observable);
Hwang et al. (2022) Section 2.5, which gives the phase-space argument for the
flat thermal law and reports numerically that it stays flat under velocity
criteria (Proposition 8.1 proves and delimits this); Hamilton (2022) on
tidally anisotropic orientations; Makarov (2025, 2026) parts I-III; Wicksell
(1925) and Groeneboom and Jongbloed (1995) for the stereological analogue.

**Outcome.** No prior statement was found of the Fourier formula for the mean
square eccentricity of the projected-angle law, of the dual-observable
uniqueness theorem for this kernel, or of the separation-binning reweighting
for this observable. The thermal selection result rests on the classical
microcanonical characterisation of the thermal law (Jeans 1919), which is
credited. The weak-continuity argument and the moment and Abel tools are
standard. One alternative route to uniqueness, through the intrinsic-angle law
of Hwang et al. (equation 16) and the projection of a rigid angle under random
orientation, was identified by a reviewer; the literatures where that
projection kernel occurs (jet-disc and pulsar spin-velocity alignments) were
not searched. This addendum is a bounded search, not a priority certificate.
