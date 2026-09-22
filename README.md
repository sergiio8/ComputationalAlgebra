# Computational Algebra

This repository is an **academic algorithms collection / coursework project** for computational algebra. It contains direct Python implementations of finite fields, polynomial rings, fast polynomial operations, Fourier transforms, and polynomial factorization algorithms. The code is intended for study and experimentation with the underlying mathematics, rather than as a packaged production library.

## Contents

### `cuerpos_finitos.py`

Foundational algebraic structures:

- `cuerpo_fp`: arithmetic in the prime field \(\mathbb{F}_p\).
- `anillo_fp_x`: univariate polynomials over \(\mathbb{F}_p\), including arithmetic, division with remainder, gcd and extended gcd, modular powers, derivatives, irreducibility tests, and helpers for tables and random polynomials.
- `cuerpo_fq`: finite extension fields \(\mathbb{F}_q\) represented as a quotient of \(\mathbb{F}_p[x]\) by a supplied modulus.
- `anillo_fq_x`: univariate polynomials over an extension field, with the analogous arithmetic and algebraic operations.

Elements are represented internally as tuples or integers according to the class that creates them. The classes expose conversion and pretty-printing methods so callers do not need to depend on those representations.

### `algoritmos_rapidos.py`

Fast and structured algorithms built on the polynomial-ring classes:

- Karatsuba-style polynomial multiplication over prime and extension fields.
- Matrix-vector products for lower, upper, and general Toeplitz matrices.
- Recursive inversion of triangular Toeplitz matrices.
- Polynomial division via Toeplitz systems.
- Cooley–Tukey FFT and inverse FFT routines over finite fields.

Importing this module attaches the fast multiplication and division routines as additional methods on the polynomial-ring classes defined in `cuerpos_finitos.py`.

### `factorizacion.py`

Polynomial factorization routines over prime and extension fields:

- Square-free factorization.
- Distinct-degree factorization.
- Equal-degree factorization.
- Multiplicity calculation.
- Cantor–Zassenhaus factorization through `fact_fpx` and `fact_fqx`.

This module depends on `cuerpos_finitos.py` and, like the coursework specification, assumes the relevant field and polynomial inputs satisfy the mathematical preconditions of each algorithm.

## Requirements and setup

The repository contains standalone Python source files and no package metadata or third-party dependency manifest. A Python installation is therefore sufficient for the code currently included; no external dependency is documented or required by the source files.

Clone the repository and run the modules from its root:

```bash
git clone https://github.com/sergiio8/ComputationalAlgebra.git
cd ComputationalAlgebra
python
```

The modules are intended to be imported by exercises or interactive sessions. For example:

```python
import cuerpos_finitos as cf

fp = cf.cuerpo_fp(5)
fpx = cf.anillo_fp_x(fp)

f = fpx.elem_de_tuple((1, 2, 1))  # 1 + 2x + x^2
g = fpx.elem_de_tuple((4, 1))     # 4 + x

product = fpx.mult(f, g)
quotient, remainder = fpx.divmod(f, g)

print(fpx.conv_a_str(product))
print(fpx.conv_a_str(quotient))
print(fpx.conv_a_str(remainder))
```

To use the faster algorithms, import the module after the foundational classes:

```python
import algoritmos_rapidos  # registers fast methods on the polynomial rings

fast_product = fpx.mult_fast(f, g)
```

## Mathematical and implementation notes

- `cuerpo_fp` validates that its modulus is prime. Extension fields require a suitable modulus polynomial supplied by the caller.
- FFT routines require a valid finite-field root of the requested power-of-two order and an input of the corresponding length.
- The factorization routines include assumptions stated in their coursework interfaces, including odd characteristic for the implemented Cantor–Zassenhaus path.
- The repository uses Spanish identifiers and comments because it originated as university coursework. The public module and class names are retained as part of that original interface.

## Scope and limitations

This is an educational implementation, not a maintained general-purpose algebra system. The repository currently has:

- no packaging configuration or command-line interface;
- no automated test suite or benchmark suite;
- no formal API stability guarantee;
- no documented support policy for Python versions;
- limited input validation beyond checks implemented by the individual algorithms.

Results should be checked against the mathematical preconditions of the selected routine before using them in other software.

## License

The project is distributed under the [MIT License](LICENSE).
