# SO(5) vector ferromagnet in a spin-2 U(1) basis

`models/contrib/so5-spin2-u1.cpp` generates a one-site infinite lattice containing
the SO(5) Hamiltonian, its broken lowerers, the embedded SO(3) generators,
the collective scalar lowerer, and finite-support convenience functions.
The only enforced symmetry is physical spin-2 `Sz:U(1)`.

Both build systems discover contributed model files automatically. From an
existing configured Autoconf build directory:

```sh
make -j4 so5-spin2-u1
./so5-spin2-u1 -o so5.lat
```

For CMake, enable `MPTK_BUILD_CONTRIB_MODELS=ON` and build the
`so5-spin2-u1` target. The model preserves the toolkit's current C++20 baseline.

## Local basis and embedding

The basis labels are the usual spin-2 values `-2`, `-1`, `0`, `1`, `2`.
The highest-weight state is the product of local state `2`.
Writing \(E_{ab}=|S^z=a\rangle\langle S^z=b|\), define

\[
Q=E_{2,2}-E_{-2,-2}=\frac{(S^z)^3-S^z}{6},\qquad
T^z=E_{1,1}-E_{-1,-1}=\frac{4S^z-(S^z)^3}{3},
\]
\[
T^+=\sqrt2(E_{1,0}+E_{0,-1}),\qquad T^-=(T^+)^\dagger,
\qquad \mathbf T^2=2P_0.
\]

Here \(P_0=1-Q^2\) projects onto the three middle levels. Thus

\[
\mathbf5\downarrow SO(2)\times SO(3)
=(+1,0)\oplus(0,1)\oplus(-1,0),\qquad S^z=2Q+T^z.
\]

The \(T\) generators refer to this reducible \(0\oplus1\oplus0\)
embedding, not the physical spin-2 SU(2).

| Local operator | Meaning | Physical Sz charge |
|---|---|---:|
| `I`, `Sp`, `Sm`, `Sz` | Standard spin-2 operators | 0, +1, -1, 0 |
| `P_plus`, `P_minus`, `P_zero` | Projectors onto \(Q=+1,-1,0\) | 0 |
| `Q` | SO(2) generator | 0 |
| `Tz`, `Tp`, `Tm` | Embedded SO(3) generators | 0, +1, -1 |
| `T2` | On-site \(\mathbf T^2=2P_0\) | 0 |
| `r_p`, `r_0`, `r_m` | Commuting broken lowerers | -1, -2, -3 |
| `r_p_dag`, `r_0_dag`, `r_m_dag` | Their adjoints | +1, +2, +3 |
| `c_minus`, `c_plus` | On-site scalar lowerer and adjoint | -4, +4 |

The three broken lowerers are

\[
r_a=E_{a,2}-(-1)^{2-a}E_{-2,-a},\qquad a=+1,0,-1.
\]

They satisfy \([r_a,r_b]=0\), \([Q,r_a]=-r_a\), and

\[
[T^z,r_a]=a r_a,\qquad
[T^\pm,r_a]=\sqrt{(1\mp a)(2\pm a)}\,r_{a\pm1}.
\]

These formulas fix the spherical-component phases used by the model.
They give the same scalar descendant ray as the Cartesian \(Z\cdot Z\)
construction, up to an overall factor and the rotation identifying the basis.

## Hamiltonian and collective operators

Let \(R_a=\sum_j r_{a,j}\). The lattice contains

| Infinite MPO | Definition |
|---|---|
| `H_SO5`, `H` | SO(5) ferromagnetic Hamiltonian |
| `Q_tot`, `Sz_tot` | Sums of the corresponding local generators |
| `Tz_tot`, `Tp_tot`, `Tm_tot` | Total embedded SO(3) generators |
| `R_p`, `R_0`, `R_m` | Collective broken lowerers |
| `T2_tot` | \((T^z_{\rm tot})^2+(T^+_{\rm tot}T^-_{\rm tot}+T^-_{\rm tot}T^+_{\rm tot})/2\) |
| `C_minus` | \(R_0^2-2R_{+1}R_{-1}\) |
| `C_plus` | \(C_-^\dagger\) |

The local bond is available as `h_SO5(j)`. With `SS(j)` denoting
physical \(\mathbf S_j\cdot\mathbf S_{j+1}\), it obeys

\[
h_{j,j+1}=1+\frac56 SS(j)-\frac19 SS(j)^2-\frac1{18}SS(j)^3.
\]

The implementation constructs the equivalent generator expression

\[
h_{ij}=-Q_iQ_j-\mathbf T_i\cdot\mathbf T_j
-\sum_{a=-1}^1(r_{a,i}^\dagger r_{a,j}+r_{a,i}r_{a,j}^\dagger).
\]

The bond spectrum is \(-1\) (multiplicity 14), \(+1\) (10), and \(+4\) (1).
This convention includes the physical bond constant: the exact energy is

\[
E_0=-(L-1)\quad\text{for an open chain},\qquad e_0=-1\quad\text{in the infinite chain}.
\]

The collective scalar contains an essential on-site contribution:

\[
C_-=\sum_i c_{-,i}
+2\sum_{i<j}\bigl(r_{0,i}r_{0,j}-r_{+1,i}r_{-1,j}-r_{-1,i}r_{+1,j}\bigr),
\qquad c_-=r_0^2-2r_{+1}r_{-1}=-3E_{-2,2}.
\]

`C_minus` lowers \(Q_{\rm tot}\) by 2 and physical \(S^z_{\rm tot}\) by 4,
while preserving the total embedded SO(3) quantum numbers. Its string terms
use an explicit total U(1) charge of -4; the scalar default of `sum_string`
would be incorrect here. `T2_tot` likewise includes both its on-site and
all-pairs contributions.

## Finite support

The generator writes an infinite lattice, as usual. For finite measurements
or MPO application, use the half-open restriction `fsup(a,b,operator)`:

```sh
mp-expectation psi 'so5.lat:fsup(0,16,H_SO5)'
mp-apply 'so5.lat:fsup(0,16,C_minus)' psi psi_lowered
```

The model also defines unit-cell functions. For example,

```text
H_finite(0){L=16}
C_minus_finite(0){L=16}
C_plus_finite(0){L=16}
T2_tot_finite(0){L=16}
Q_tot_finite(0){L=16}
H_select_finite(0){L=16,lambda=1}
```

All collective operators in the table have a `<name>_finite(cell){L}` wrapper;
`H_finite` is an additional alias for `H_SO5_finite`.
The function's cell argument translates its support to `[cell,cell+L)`.
Use integer \(L\geq1\). `H_finite` contains the \(L-1\) internal bonds;
`T2_tot_finite` is the square of the total \(T\) **inside this interval**.
Collective operators are represented by infinite triangular MPOs for finite
restriction. Their formal infinite sums need not have finite thermodynamic
expectations. Use `H_SO5` for the original infinite-system Hamiltonian.

## Explicit scalar-state construction

For a four-site example:

```sh
so5-spin2-u1 -o so5.lat
mp-construct -l so5.lat -o omega.psi --finite -r 4 2
mp-apply 'so5.lat:C_minus_finite(0){L=4}' omega.psi scalar.psi
mp-normalize scalar.psi
mp-expectation --real scalar.psi 'so5.lat:H_finite(0){L=4}'
mp-expectation --real scalar.psi 'so5.lat:T2_tot_finite(0){L=4}'
mp-expectation --real scalar.psi 'so5.lat:Q_tot_finite(0){L=4}'
mp-info --entropy --base2 scalar.psi
```

The expected results are \(E=-3\), \(T_{\rm tot}^2=0\), \(Q_{\rm tot}=2\),
and half-chain entropy \(2.25105939349207\) bits. Before normalization,

\[
\|C_-|\Omega_L\rangle\|^2=3L(2L+1).
\]

More generally, apply and normalize after each of \(L/4\) applications for

\[
L=4t,\qquad C_-^t|\Omega_L\rangle
\propto|L;q=L/2,\ell=0,m=0\rangle.
\]

Repeated exact MPO application can grow the MPS bond dimension. For large
systems, use controlled compression or the variational route below.

## Selection by minimization

The infinite-MPO function `H_select{lambda=1}` is

\[
H_{\rm select}=H_{SO(5)}+\lambda\,\mathbf T_{\rm tot}^{\,2}.
\]

For finite DMRG, pass this triangular MPO directly. The finite DMRG driver
supplies the boundary contraction; do not wrap its `-H` argument in `fsup`.
For \(L=4t\), fix physical \(S^z_{\rm tot}=L\):

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mp-construct -l so5.lat -o trial.psi --finite -r 8 2:0
mp-dmrg -w trial.psi -H 'so5.lat:H_select{lambda=1}' -m 128 -s 32 \
  --pre rsvd --pre-factor 0.1 --post rsvd --post-factor 0.1 \
  --miniter 2 --maxiter 8
mp-expectation --real trial.psi 'so5.lat:H_finite(0){L=16}'
mp-expectation --real trial.psi 'so5.lat:T2_tot_finite(0){L=16}'
mp-info --entropy --base2 trial.psi
```

Repeat batches of half-sweeps and check convergence of the environments and
observables. Prefer short local Lanczos solves and more sweeps; the example
uses the adaptive default local tolerance rather than a uniformly stringent
local tolerance. The example uses `rsvd` post-expansion; `mixing` is also available.
Keep OpenBLAS at one thread for these calculations.

In the local 16-site validation, the settings above, with normalization after
each batch, converged after seven batches (224 half-sweeps). The final entropy
was 4.61465522888757 bits, its error against the exact scalar spectrum was
7.1e-14 bits, and the overlap with explicit lowering rounded to one. Convergence
was confirmed at two consecutive measurement points; a fixed sweep count alone
is not a convergence criterion. Initial runs with `fullsvd` post-expansion and
50--100 local Lanczos iterations are retained in the calculation record as
exploratory results, not as recommended settings.

The stopping thresholds in this validation were intentionally stringent.
For routine work, consider the measured truncation error and whether further
sweeps or larger bond dimensions meaningfully change the plotted observables.
Already at 96 half-sweeps, the entropy error was only -1.50e-6 bits, far below
the finite-size entropy figure's graphical resolution.

The scalar ground state is unique in this magnetization sector when
\(\lambda>0\): \(T_{\rm tot}^2=0\) implies \(m=0\), and
\(S^z_{\rm tot}=2q+m=L\) implies \(q=L/2\). The SO(5) ground irrep
contains this \((q,\ell)\) exactly once. The selector commutes with the
original Hamiltonian and annihilates the target, so no zero-penalty
extrapolation is needed. It changes the excited-state spectrum and is meant
for selecting finite ground-state vectors, not for calculating the original
model's Goldstone dispersion.

Converge bond dimension and sweeps, check the original energy and
\(T_{\rm tot}^2\), and compare the Schmidt spectrum with the analytic one.
Energy convergence alone is insufficient for the unselected degenerate model.

## Verification

The contributed integration suite is `tests/suites/contrib/so5-spin2-u1.yaml`:

```sh
scripts/mptk-test --bin-dir /path/to/build contrib/so5-spin2-u1
```

It checks highest-weight energies, the on-site part of \(C_-\), equivalence
of `fsup` and the convenience functions, scalar-state quantum numbers and
entropy, and minimization against explicit lowering. The spin-2 Hamiltonian
is Appendix B of [arXiv:2503.22468v1](https://arxiv.org/html/2503.22468v1).

The generator was compiled with the Autoconf build against public MPToolkit
commit `7c850713fc057e1a75009b661de178c75de1e9c9`. All seven integration
checks passed, including four-site DMRG versus explicit lowering and
eight-site repeated lowering with half-chain entropy
\(3.324362673901629\) bits. These are small-system checks; large-system
convergence and infinite DMRG have not been tested here.
