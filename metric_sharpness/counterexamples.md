# A certified metric boundary for prefix-cone thermal ordering

Research note, 26 October 2026. All indices below are one-based. This note proves a boundary of the abstract reversible cone-preserving theorem. It makes no claim that the Rényi cutoff at order two is sharp.

## 1. Result and scope

**Proposition.** There exists a five-state continuous-time Markov generator with all transition rates strictly positive, detailed balance with a strictly ordered Gibbs energy spectrum, and preservation of the prefix-mean cone
\[
\mathcal C=\left\{b:b_j\ge\frac1{j-1}\sum_{i<j}b_i,
\quad j=2,\ldots,N\right\},
\]
for which two thermal preparations on the heating side of the bath have a total-variation crossing. Consequently these hypotheses do not guarantee conventional heating order for every convex \(f\)-divergence.

The example is a single time-independent reversible generator, so it also belongs to the finite-protocol class. It does not contradict the proved heating order for convex \(f\) with convex derivative: total variation is outside that sufficient class. The same example has conventional chi-square order by that theorem. Its two evolved states are therefore incomparable in thermomajorization at the exhibited time.

This result concerns **heating**. It does not establish a cooling TV crossing, a cooling Rényi crossing above order two, or a heating Rényi crossing below order two. Those questions remain unresolved by this branch. Nor is the generator asserted to arise from equal pair couplings and a concave bosonic jump kernel.

## 2. Four-state mechanism

For a row generator acting on observables, set
\[
G_4=\begin{pmatrix}
-113/552&203/1104&1/48&0\\
145/184&-313/368&0&1/16\\
5/8&0&-11/16&1/16\\
0&7/8&1/8&-1
\end{pmatrix}. \tag{1}
\]
It has stationary law \(\bar\pi=(60,14,2,1)/77\). The detailed-balance equalities follow directly, for example
\(60(203/1104)=14(145/184)\), and similarly on the other three edges. The graph is connected, so the chain is irreducible despite its two absent edges.

Let \(Q_{ij}=G_{ij}-G_{1j}\) for \(i,j\ge2\). If \(h_j=b_j-b_1\), the lower-triangular matrix \(T\) maps \(h\) to
\[
d_j=h_j-\frac1{j-1}\sum_{i=2}^{j-1}h_i.
\]
The cone condition is \(d\ge0\). Direct exact multiplication gives
\[
H_4=TQ T^{-1}=\begin{pmatrix}
-373/368&0&1/16\\
0&-11/16&1/32\\
3/4&0&-25/24
\end{pmatrix}. \tag{2}
\]
This is Metzler. Thus \(e^{tG_4}\mathcal C\subseteq\mathcal C\) for all \(t\ge0\).

Use bath inverse temperature \(\beta_b=1\) and energies
\[
E_i=\log(60/w_i),\qquad (w_1,w_2,w_3,w_4)=(60,14,2,1).
\]
At integer preparation inverse temperatures the Gibbs probabilities are rational:
\[
p_{\beta,i}=\frac{w_i^\beta}{\sum_jw_j^\beta}.
\]
Take the warm and cold preparations \(\beta_w=3\), \(\beta_c=4\). Both heat towards the bath. Define
\[
k=\frac{9992521}{10^7},\qquad
f_k(x)=(x-k)_+-(1-k)_+,
\quad D_k(p\Vert\bar\pi)=\sum_i\bar\pi_i f_k(p_i/\bar\pi_i).
\tag{3}
\]
This is a convex generator with \(f_k(1)=0\). The exact certificate described below proves
\[
D_k(p_4\Vert\bar\pi)-D_k(p_3\Vert\bar\pi)
>0.00962822924,
\tag{4}
\]
whereas at \(t_*=1047/100\), with row-distribution evolution \(p_\beta(t)=p_\beta e^{tG_4}\),
\[
6.76704\times10^{-8}
<D_k(p_3(t_*)\Vert\bar\pi)-D_k(p_4(t_*)\Vert\bar\pi)
<6.76705\times10^{-8}. \tag{5}
\]

The likelihood ratios of states \(1,2,4\) exceed \(k\) for both evolved preparations, while that of state \(3\) lies below \(k\). Thus the strict sign in (5) is simply the sign of the difference of the probabilities of the fixed event \(\{1,2,4\}\); no uncertain sorting decision enters the certificate.

## 3. Five-state total-variation example

Adjoin the highest-energy weight
\[
w_5=\frac{5763}{100000},\qquad
w=(60,14,2,1,5763/100000),\qquad
\pi_i=\frac{w_i}{\sum_jw_j}.
\tag{6}
\]
Equivalently,
\[
\pi=\frac1{7705763}(6000000,1400000,200000,100000,5763).
\]
Use the same energy prescription \(E_i=\log(60/w_i)\), bath \(\beta_b=1\), and preparations \(\beta_w=3,\beta_c=4\). The energies are strictly increasing and all probabilities have full support.

The generator is the following finite rational prescription:
\[
G_5=\begin{pmatrix}G_4&0\\0&0\end{pmatrix}
+\frac1{100}(\mathbf1\pi^T-I_5). \tag{7}
\]
The block-diagonal part is reversible with respect to \(\pi\), since its active block has the same stationary ratios as \(\bar\pi\). The added Gibbs-reset generator is also reversible. Their sum is therefore reversible, and every off-diagonal rate obeys
\[
(G_5)_{ij}\ge\frac{\pi_j}{100}
\ge\frac{5763}{770576300}>0. \tag{8}
\]
This proves irreducibility and strict positivity without a limiting argument.

### 3.1 Exact cone-preservation check

For the block-diagonal generator before the reset, the five-state transformed quotient is
\[
H_{5,0}=\begin{pmatrix}
-373/368&0&1/16&0\\
0&-11/16&1/32&0\\
3/4&0&-25/24&0\\
1/8&5/24&7/32&0
\end{pmatrix}. \tag{9}
\]
The new last row is the explicit check that the spectator embedding preserves the prefix-mean cone; that property is not presumed for an arbitrary direct sum.

The reset generator has quotient \(-I\), because for \(i,j\ge2\),
\[
(\pi_j-\delta_{ij})-\pi_j=-\delta_{ij}.
\]
It therefore contributes \(-I/100\) after conjugation by \(T\). Consequently
\[
H_5=H_{5,0}-I_4/100
=\begin{pmatrix}
-9417/9200&0&1/16&0\\
0&-279/400&1/32&0\\
3/4&0&-631/600&0\\
1/8&5/24&7/32&-1/100
\end{pmatrix}, \tag{10}
\]
which is Metzler. This proves the required cone invariance exactly.

### 3.2 Certified TV crossing

With \(\operatorname{TV}(p,\pi)=\frac12\sum_i|p_i-\pi_i|
=\sum_i(p_i-\pi_i)_+\), the rational initial states give
\[
\operatorname{TV}(p_4,\pi)-\operatorname{TV}(p_3,\pi)
>0.00962823010. \tag{11}
\]
At \(t_*=1047/100\), exact rational uniformisation bounds give
\[
6.01768\times10^{-8}
<\operatorname{TV}(p_3(t_*),\pi)
-\operatorname{TV}(p_4(t_*),\pi)
<6.01769\times10^{-8}. \tag{12}
\]
Hence continuity guarantees a finite crossing before \(t_*\).

For both preparations the positive event in the TV expression at \(t_*\) is exactly \(\{1,2,4\}\), and the negative event is \(\{3,5\}\). The two most delicate ratio inequalities are those of the cold preparation:
\[
0.99999991312647468<r_{4,3}(t_*)<0.99999991312647470<1,
\]
\[
1<1.00000005348966451<r_{4,4}(t_*)<1.00000005348966453.
\]
They are separated from one by more than \(5\times10^{-8}\), far above the certified arithmetic uncertainty.

The spectator changes the reference normalisation so that a small interval of four-state likelihood ratios near \(k\) straddles the five-state TV threshold at one. The reset supplies strictly positive rates while retaining the relevant structure. More precisely, if \(B=\operatorname{diag}(G_4,0)\) and \(\Pi=\mathbf1\pi^T\), then \(B\Pi=\Pi B=0\), and
\[
e^{t[B+\eta(\Pi-I)]}=\Pi+e^{-\eta t}(e^{tB}-\Pi).
\]
Thus the reset multiplies every deviation from \(\pi\) by a common positive factor. It preserves any TV ordering of this embedded model. The certificate nevertheless checks the full generator in (7) directly.

## 4. Fully explicit error bounds

The verification program `certify_counterexamples.py` uses Python's exact `Fraction` arithmetic for every admissibility and sign decision. Decimal values in `exact_counterexamples.json` are solely readable renderings.

For either generator, choose the rational uniformisation rate \(\lambda=2\) and
\[
P=I+G/2,\qquad x=2t_*=1047/50=20.94.
\]
The program checks exactly that \(P\) is stochastic and entrywise nonnegative. For each rational initial row \(p\), it computes
\[
S_{100}=\sum_{n=0}^{100}\frac{x^n}{n!}pP^n,
\qquad A_{100}=\sum_{n=0}^{100}\frac{x^n}{n!}.
\]
Because \(P^n\) is stochastic, each component of the omitted vector series is between zero and the scalar remainder bound
\[
R=\frac{x^{101}}{101!}
\frac1{1-x/102}
<3.501\times10^{-27}. \tag{13}
\]
The inequality follows by bounding all subsequent term ratios by \(x/102<1\). The same bound encloses the omitted terms of \(e^x\). Thus, componentwise,
\[
\frac{S_{100,i}}{A_{100}+R}
\le (p e^{t_*G})_i
\le\frac{S_{100,i}+R}{A_{100}}. \tag{14}
\]
Every endpoint is rational. The largest resulting probability-interval width for either model is less than \(5.015\times10^{-36}\). The program first proves that the relevant threshold sets are common to the two preparations, then sums the lower bounds for the warm state and upper bounds for the cold state on that set. The resulting lower bounds are strictly positive and give (5) and (12). The opposite endpoint choices give the displayed upper bounds.

To reproduce all certificates:

```sh
python physics_research/four_star_research_2026-10-26/metric_sharpness/certify_counterexamples.py
```

The script also proves detailed balance, row sums, all off-diagonal signs, the transformed cone inequalities, the strictly ordered rational weights, and the initial ordering. There is no eigenvalue computation, numerical integration, floating-point comparison, or unbounded series truncation in the certificate.

## 5. Interpretation and limitations

The metric boundary is substantive: preservation of the energy prefix cone and detailed balance do not extend the heating result to all convex \(f\), and they do not give a TV theorem covering both cooling and heating. Conventional energy-distribution order can coexist with a transient TV inversion. This is why a theorem stated for the proved derivative classes must not be promoted to thermomajorization order.

The example has modest finite parameters, but its inversion is small. It is an exact counterexample and a scope statement, not evidence of a practically useful speedup. It supplies neither a strong-Mpemba cancellation nor a universal quantitative acceleration result. The particular four-state construction and spectator normalisation are tuned to expose the sign failure.

Metric dependence and thermomajorization are established concepts. The object requiring novelty assessment is specifically a thermal heating crossing **inside the reversible prefix-cone preserving class**, with the explicit positivity and cone certificates above. The bounded source comparison is summarised in `PRIOR_ART.md` at the package root.

## 6. Late-time return and at least two crossings

The five-state inversion is transient. The following analytical argument supplements the rational finite-time certificate; it does not require another numerical calculation.

Put \(B=\operatorname{diag}(G_4,0)\), \(\Pi=\mathbf1\pi^T\), and \(\eta=1/100\). As already noted, \(B\Pi=\Pi B=0\), so
\[
e^{tG_5}=\Pi+e^{-\eta t}(e^{tB}-\Pi). \tag{15}
\]
The Gershgorin right edges of the three rows of \(H_4\) in (2) are, respectively,
\[
-\frac{175}{184},\qquad -\frac{21}{32},\qquad-\frac7{24}.
\]
Every nonzero eigenvalue of \(G_4\) therefore has real part at most \(-7/24\). Reversibility makes the generator diagonalisable, so its active-block relaxation is \(O(e^{-7t/24})\), with no polynomial Jordan factor.

Let \(q_\beta=p_{\beta,5}\) denote the initial probability of the spectator, and let \(\bar\pi=(60,14,2,1)/77\). Since \(B\) does not transfer mass between the two blocks,
\[
p_\beta e^{tB}
=((1-q_\beta)\bar\pi,q_\beta)+O(e^{-7t/24}). \tag{16}
\]
Also \(\pi=((1-\pi_5)\bar\pi,\pi_5)\), so (15) implies
\[
p_\beta e^{tG_5}-\pi
=(\pi_5-q_\beta)e^{-\eta t}(\bar\pi,-1)
+O(e^{-(\eta+7/24)t}). \tag{17}
\]
The displayed signed vector has total variation one. Consequently, for every fixed \(\beta>1\),
\[
\operatorname{TV}(p_\beta e^{tG_5},\pi)
=(\pi_5-q_\beta)e^{-\eta t}
+O(e^{-(\eta+7/24)t}). \tag{18}
\]
The highest-level Gibbs probability is strictly decreasing in inverse temperature:
\[
q_\beta'=-q_\beta\bigl(E_5-\langle E\rangle_{p_\beta}\bigr)<0,
\]
because \(E_5\) is strictly larger than the other energies and all probabilities are positive. Thus
\(q_4<q_3<q_1=\pi_5\), and subtraction in (18) gives
\[
\operatorname{TV}(p_4 e^{tG_5},\pi)
-\operatorname{TV}(p_3 e^{tG_5},\pi)
=(q_3-q_4)e^{-\eta t}
+O(e^{-(\eta+7/24)t})>0
\]
for all sufficiently large finite \(t\).

The same difference is positive initially by (11) and negative at \(t_*\) by (12). Continuity therefore supplies **at least two distinct sign-changing crossings**: one before \(t_*\) and one afterwards. This does not assert that there are exactly two crossings or locate the second one numerically.

The block-mass mode in (17) is the unique slow mode: \(B\) has two stationary block modes and three modes bounded above by \(-7/24\), while reset leaves one stationary mode and moves the other to \(-\eta\). Its signed thermal coefficient \(\pi_5-q_\beta\) is strictly increasing in \(\beta\) and vanishes only at the bath. Its absolute value therefore has conventional ordering on either side of the bath. The generator exhibits a transient TV inversion without a slow-mode inversion or a nontrivial strong-Mpemba cancellation.

This late-time derivation is also present in the coordinating manuscript's reset identity and TV asymptotic equations and has been independently checked by its reviewer. The certificate files remain unchanged by this analytical addition.
