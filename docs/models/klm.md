# Kondo lattices with arbitrary local spin and singlet pairing

`klm-u1su2` conserves electron number and total spin. `klm-su2` conserves total
spin and supports onsite singlet pairing. Both represent one conduction-electron
orbital and one localized moment per site, repeated along a chain; they are
Kondo lattices, not single-impurity models.

From a configured build directory:

```sh
make klm-u1su2 klm-su2
./klm-u1su2 --Spin 1 -o klm-number.lat
./klm-su2 --Spin 1 -o klm-paired.lat
```

CMake discovers both generators automatically. `--Spin` (or `-S`) accepts
nonnegative integers and half-integers written as decimals, e.g. `0.5`, `1`,
`1.5`. The default remains `0.5`. Like the other spin generators, these use the
existing `half_int` stream conversion, which rounds other decimal inputs to
the nearest half-integer (for example, `0.7` becomes `0.5`). The site constructor
rejects negative spins. There is no Kondo-specific conversion tolerance or spin
cap. `klm-cylinder-u1su2` also accepts this option. The separate `klm-u1u1`
generator is unchanged.

## Local basis and compatibility

The physical site dimension is `4*(2*S+1)`. The reduced basis is:

| Label | Electron number | Total site spin |
|---|---:|---:|
| `empty` | 0 | S |
| `double` | 2 | S |
| `lower` | 1 | S - 1/2 |
| `upper` | 1 | S + 1/2 |

At S=1/2, the original labels `singlet` and `triplet` replace `lower` and
`upper`, with the original ordering, phases, and operator normalizations.
At S=0, `lower` is absent and the site reduces to a conduction orbital with
no local moment. The labels still distinguish empty and double occupancy in
the SU(2)-only model, even though those multiplets have identical quantum numbers.

`Sc` is the conduction spin, `Sf` the localized moment, and `S = Sc + Sf` the
total spin. In particular, `S` does not mean the local moment alone. `N_S` is
the projector onto single electron occupancy, not the number or magnitude of
localized spins. `N_H = I - N_S`. Empty and doubly occupied sites still carry
the localized spin.

The conduction spin uses the normalization `Sc = (1/2) c^dagger sigma c`,
where `sigma` denotes the Pauli matrices. If another Hamiltonian instead
defines its conduction spin as `c^dagger sigma c`, its Kondo coupling must
be multiplied by two when used as the coefficient of `H_K` or `ScSf` here.

The electron spinors `C`/`CH` and the spin vectors are recoupled using the
tensor-product machinery. `ScSf` has eigenvalues `-(S+1)/2` in the lower
multiplet, `S/2` in the upper multiplet, and zero at empty/double occupancy.
Thus at S=1 these are -1 and +1/2, rather than -3/4 and +1/4.

The C++ constructors are `KondoU1SU2(Spin, charge_name, spin_name)` and
`KondoSU2(Spin, spin_name)`. The old `KondoU1SU2()` and
`KondoU1SU2(charge_name, spin_name)` calls remain valid with S=1/2.

## Hamiltonian terms

The chain generators provide:

| Operator | Definition |
|---|---|
| `H_t` | minus nearest-neighbor electron hopping plus its Hermitian conjugate |
| `H_t2` | minus next-nearest-neighbor hopping plus its Hermitian conjugate |
| `H_K` | sum of `Sc . Sf` on each site |
| `H_J` | sum of `Sf(i) . Sf(i+1)` |
| `H_U` | sum of `n_up*n_down` |
| `H_N` | total electron number |
| `H_Delta` (SU2 only) | sum of `PairH + Pair` |
| `H_DeltaI` (SU2 only) | sum of `i*(PairH - Pair)` |

Here `PairH = c_up^dagger*c_down^dagger` and `Pair = c_down*c_up`.
For Delta = DeltaR + i*DeltaI, the Hamiltonian expression is

```text
t*H_t + JK*H_K + J*H_J + U*H_U - mu*H_N
  + DeltaR*H_Delta + DeltaI*H_DeltaI
```

Substitute numerical couplings when passing an expression to a solver. For
example, a two-site total-spin singlet with spin-1 local moments can be run as:

```sh
./mp-construct -l klm-paired.lat -o psi --finite '(0)empty(1)empty(0)'
./mp-dmrg -w psi -H 'klm-paired.lat:H_t+H_K+0.2*H_J-0.3*H_N+0.4*H_Delta' -m 40x8
```

The chemical-potential convention above is **-mu N**; use the opposite sign
if defining the Hamiltonian with +mu N. The lattice file exposes terms; it
does not select coupling constants or a default Hamiltonian.

The legacy local eta-pairing operators are preserved:

```text
Ep = 2*PairH
Em = -2*Pair = -adjoint(Ep)
```

Consequently the normalized real pairing term is `(Ep-Em)/2`, not `Ep+Em`.
The new `Pair` and `PairH` names avoid this normalization/sign trap. They are
also available in the U(1)xSU(2) site as charged operators for correlations,
but a Hamiltonian with a nonzero pairing source requires the SU(2)-only model.

### Signs of the couplings

For the convention above, positive `t` multiplies `H_t`: the minus sign in
`-t sum(c_i^dagger c_{i+1} + h.c.)` is already included in `H_t`. Positive
`JK` and `J` are antiferromagnetic. For example, at local spin 1, positive
`JK*ScSf` lowers the onsite total-spin-1/2 multiplet (eigenvalue `-JK`)
relative to the total-spin-3/2 multiplet (eigenvalue `JK/2`). The two-site
local-moment singlet has `H_J = -2`.

`H_Delta` has no built-in minus sign. A term
`+Delta*(c_up^dagger*c_down^dagger + c_down*c_up)` therefore requires
`+Delta*H_Delta`; `-Delta*H_Delta` describes the opposite pairing phase.
Reversing the order of either pair of fermion operators also reverses its
sign. For a single uniform pairing source, changing its overall sign is a
global gauge transformation (`c -> i*c`): energies are unchanged, but
anomalous expectations such as `Pair` change sign. Relative pairing phases
must still be kept consistent when comparing different regions or sources.

### Different Kondo couplings at the ends of a finite chain

The `mp-dmrg -H` interface takes an infinite triangular MPO and contracts it
with the finite wavefunction's open boundaries. To use bulk coupling `JK`
and end couplings `JL` and `JR` on sites `0,...,N-1`, add the differences
from the bulk coupling:

```text
Hbulk + sum_unit(sites=N, (JL-JK)*ScSf(0) + (JR-JK)*ScSf(N-1))
```

Here `Hbulk` includes `JK*H_K`. These names are placeholders: substitute
numerical values or expand shell variables before passing the expression
to MPToolkit. The subtraction is essential because `H_K` already includes
the Kondo interaction at both ends.

`sum_unit(sites=N, op)` repeats `op` every **N sites** on the infinite
lattice. Restricting this construction to the finite chain `[0,N)` leaves
the two intended end terms. It does not impose periodic boundary conditions
on the finite wavefunction. Match `N` to the chain length; using the same
operator on a longer chain repeats the modifications in its interior.

Both chain generators have a one-site lattice unit cell, so `ScSf(11)`
means site 11, the right end of a 12-site chain. More generally,
`Op(cell)[site]` selects a site within a lattice unit cell; the cell index
still refers to the generator's unit cell, even inside `sum_unit(sites=N,...)`.
The `sites=` argument must precede the operator expression.
`sum_unit(ScSf(0),12)` is not the supported syntax. A bare `ScSf(0)` is a
finite-support operator, so it cannot be added directly to `H_K` in `-H`.

For example, start with a 12-site total-spin singlet and use `t=0.1`, `J=1`,
`JK=0.325`, `Delta=+0.2`, and `mu=U=0`, using the `Sc` normalization above.
For equal end couplings `Jedge=0.315`:

```sh
./klm-su2 --Spin 1 -o klm-paired.lat
./mp-construct -l klm-paired.lat -o psi_s0 --finite -r 6 '(0)empty(1)empty(0)'
Jedge=0.315
H12="0.1*H_t+H_J+0.325*H_K+0.2*H_Delta+sum_unit(sites=12,(${Jedge}-0.325)*(ScSf(0)+ScSf(11)))"
./mp-dmrg -w psi_s0 -H "klm-paired.lat:${H12}" -m 100x20
```

Choose the state count and sweep schedule by convergence checks. Changing
`Jedge` from `0.310` to `0.325` changes the correction from `-0.015` to zero.
Reassign `H12` after changing `Jedge` so the shell expands the new value.
For unequal end couplings, use the two separate differences in the general
expression above. If the desired pairing term is negative, replace
`+0.2*H_Delta` with `-0.2*H_Delta` explicitly.

For measurements, `mp-expectation` takes a finite-support operator. Use
`fsup` to restrict the same infinite Hamiltonian, or give local operators
directly:

```sh
./mp-expectation psi_s0 "klm-paired.lat:fsup(0,12,${H12})"
./mp-expectation psi_s0 'klm-paired.lat:ScSf(0)'
./mp-expectation psi_s0 'klm-paired.lat:ScSf(11)'
```

`fsup(a,b,operator)` uses a half-open interval of **lattice unit cells**:
the first is included and the second excluded. For these one-site cells,
`fsup(0,12,...)` keeps the terms contained in sites 0 through 11, with no
bond joining site 11 to site 0. It returns a finite-support MPO; use it for
measurement, not as the outer expression passed to `mp-dmrg -H`.

## Why no separate Z2 label?

For L sites of fixed local spin S, the central element of total-spin SU(2)
gives

```text
(-1)^(2*J_total) = (-1)^N_e * (-1)^(2*S*L).
```

Thus electron parity is already fixed by the total-spin sector and the known
local-moment offset. For integer S it is exactly the integer/half-integer
distinction; for half-integer S there is an additional length-dependent sign.
A separate Z2 quantum number would be redundant here. The fermion-parity
operator `P` and fermionic commutation metadata remain essential and are
retained in `klm-su2` for the Jordan-Wigner signs.

These generators are uniform lattices. A finite magnetic segment attached
to electron-only superconducting leads requires a mixed-site construction
(using S=0 on the leads), not just a different coupling in this uniform
generator. Spatially varying couplings can be assembled from the local
operators, but no lead geometry is generated by these commands.
