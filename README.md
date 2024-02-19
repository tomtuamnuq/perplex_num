# perplex_rs
## Overview
`perplex_rs` is a Rust crate that provides an implementation of perplex numbers, based on the [nalgebra](https://docs.rs/nalgebra/latest/nalgebra/) crate. This library allows for the representation and manipulation of perplex numbers, offering both a struct-based approach and a matrix form representation.

## Perplex Numbers
Perplex numbers, also known as **split-complex**, **double** or **hyperbolic** numbers, are an extension of the real numbers $\mathbb R$ by introducing a new element $j$ with the property that $j^2=1$, which is distinct from the imaginary unit $i$ in complex numbers ($i^2=-1$). $j$ is referred to as **hyperbolic unit**. A perplex number $z$ is expressed as:
$$z=t + x j \quad t,x \in \mathbb R \quad j^2=1$$
The perplex numbers have applications in various fields for instance special relativity. In the context of Minkowski space, which is used in special relativity, the variables $t$ and $x$ typically represent **time** (the real part of $z$) and **space** (a spatial coordinate - the hyperbolic part of $z$). A thorough description of hyperbolic numbers in this regard can be found in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6). They form a two-dimensional commutative algebra over the real numbers, similar to the complex plane, but with a different geometric interpretation due to the hyperbolic unit: 

**TODO** add an image of perplex number plane made with Plotters

### Operations
The perplex numbers form an algebraic ring with addition and multiplication (see [Wikipedia](https://wikipedia.org/wiki/Split-complex_number) for a definition in terms of abstract algebra). Let $z_1=t_1+x_1j$ and $z_2=t_2+x_2j$ be two perplex numbers:

- **Addition** is component-wise: $$z_1+z_2 = (t_1+t_2) + (x_1+x_2)j$$
- **Multiplication** is given by: $$z_1 z_2 = (t_1t_2 + x_1x_2 ) + (t_2x_1 + t_1x_2)$$

Let $z=t+xj$ be a perplex number:
- **Conjugate** of a perplex number is: $$\bar z = t - xj$$
- **Modulus** (or norm) is defined by (the square root of) the quadratic form $N$: $$N(z)=t^2-x^2$$

**TODO** space time and light-like

Perplex numbers include elements called null vectors or zero divisors, which are of the form $x+xj$ or $x-xj$, with $x\not=0$. These are exactly the non-zero elements with a modulus of zero. These numbers are not invertible, meaning they do not have a multiplicative inverse within the set of perplex numbers.

### Matrix Representation


## Features
- `Perplex` struct implements all common mathematical operations via `std::ops` traits
- 

## Usage

## Installation
`cargo add perplex_rs` or add this to `Cargo.toml`:

```toml
[dependencies]
perplex_rs = "0.1"
```

## Bibliography
- [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6)
- [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643)
- [Introduction to Hybrid Numbers](https://doi.org/10.1007/s00006-018-0833-3)

