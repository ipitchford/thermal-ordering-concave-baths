# Exact nonconcave separator with distinct levels

This is an independent audit example. It tests the all-dimension quantifier in the universal prefix-cone classification. It does **not** demonstrate an actual Mpemba crossing.

Let J(d)=1+d² on [0,1]. For any a,b>0 with a+b≤1,

\[
J(a)+J(b)-J(a+b)=1-2ab\ge\frac12.
\]

Hence every zero-bath triplet has a positive lower off-diagonal prefix-coordinate entry. Nonconcavity cannot universally be detected by a zero-bath triplet test.

Take N=16, h=1/4, and the distinct spectrum

\[
E=(0,\tfrac14,\tfrac12-\tfrac6{1000},\tfrac12-\tfrac5{1000},\ldots,
\tfrac12+\tfrac6{1000},\tfrac34).
\]

There are thirteen middle levels. The fixed-length increment is affine:

\[
A(z)=J(z+\tfrac14)-J(z)=\frac z2+\frac1{16}.
\]

Their arguments z=E_i−h have mean exactly x=1/4, so their summed A-values equal thirteen times A(1/4) exactly. Thus the zero-bath transformed entry agrees exactly with the coalesced-limit expression, despite the strictly distinct levels:

\[
H^A_{16,2}
=\frac{2J(1/4)+13A(1/4)-15A(1/2)}{30}
=-\frac1{240}.
\]

Now fix the finite bath β_b=16 log2 and h_b(d)=J(d)/(exp(β_b d)−1). Symmetry of the thermal generator gives

\[
H^S_{N,2}=\frac{N}{2(N-1)}
\left[h_b(E_N-E_2)-h_b(E_N-E_1)\right].
\]

Every middle thermal contribution cancels from this specific entry. Consequently the only surviving exponentials are exp(β_b/2)=2⁸ and exp(3β_b/4)=2¹²:

\[
H^S_{16,2}
=\frac8{15}\left[\frac{5/4}{255}-\frac{25/16}{4095}\right]
=\frac{1007}{417690}.
\]

The full entry is therefore the exact rational number

\[
\boxed{H_{16,2}=-\frac1{240}+\frac{1007}{417690}
=-\frac{5867}{3341520}<0.}
\]

Most other rates of this finite-bath model are irrational, because the middle spacings involve thousandths. No claim that the whole generator is rational is made. The exact rational entry follows from the two cancellations above. Direct assembly of the full floating-point generator gives −0.0017557877851876458, consistent with the exact result.

One can strengthen the triplet comparison at this same finite bath. On [0,1],

\[
\frac{d}{dd}\log(1+d^2)=\frac{2d}{1+d^2}\le1,
\qquad
\frac{d}{dd}\log(e^{\beta_b d}-1)
=\frac{\beta_b}{1-e^{-\beta_b d}}>\beta_b>1.
\]

Thus h_b is strictly decreasing. Every thermal triplet off-diagonal prefix entry is positive, and its downward lower entry is at least 1/4. All triplets preserve the cone at this bath, while the full 16-state generator does not. This final comparison concerns cone preservation throughout; no thermal distance reversal is inferred.

## Smaller independently checked example supplied by the physical branch

After the 16-state example was recorded, the physical branch supplied the smaller distinct spectrum

\[
E=(0,1/2,501/1000,502/1000,503/1000,1).
\]

For h=1/2, the increment becomes A(z)=z+1/4. Direct transformation, independently implemented in `independent_checks.py`, gives

\[
H^A_{6,2}
=\frac{2(5/4)+(251/1000+252/1000+253/1000)-5(3/4)}{10}
=-\frac{247}{5000}.
\]

At β_b=16 log2, the same endpoint cancellation yields

\[
H^S_{6,2}=\frac35\left[\frac{5/4}{255}-\frac2{65535}\right]
=\frac{1277}{436900},
\]

and therefore

\[
\boxed{H_{6,2}=-\frac{1015293}{21845000}<0.}
\]

Every width-at-most-one triplet still has the positive downward margin and decreasing thermal kernel established above. This six-state example is easier to display in a manuscript. It is not claimed to be minimal. Its source is the physical branch; this review independently checked its arithmetic and finite-bath extension.
