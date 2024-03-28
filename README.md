# perplex_num
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE-MIT) [![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](./LICENSE-APACHE) [![minimum rustc 1.76](https://img.shields.io/badge/rustc-1.76+-red.svg)](https://rust-lang.github.io/rfcs/2495-min-rust-version.html)[![codecov](https://codecov.io/gh/tomtuamnuq/perplex_num/graph/badge.svg?token=EGEEP9PBHX)](https://codecov.io/gh/tomtuamnuq/perplex_num)

## Overview
`perplex_num` is a Rust crate that provides an implementation of perplex numbers, based on the numerical abstractions of the [num_traits](https://docs.rs/num-traits) crate. This library supports various mathematical functions such as `pow`, `sqrt`, `exp`, `ln`, `sinh`, `sin`, `cosh`, and `tan`. Additionally, the crate offers a hyperbolic polar form for representing and manipulating numbers in the hyperbolic plane, as well as a matrix form representation feature based on [nalgebra](https://docs.rs/nalgebra).

For an **in-depth explanation** (including visualizations) of perplex numbers and how they integrate with the crate's modules, see the [Perplex Number Description](https://github.com/tomtuamnuq/perplex_num/blob/main/Perplex.md) in the repository.

## Features
- The `Perplex` struct is equipped with a comprehensive set of common mathematical operations, courtesy of `std::ops` and `num_traits`.
- Emulating the functionality of `nalgebra::Complex`, the `Perplex` struct mirrors most functions found in the [num_complex](https://github.com/rust-num/num-complex) crate, maintaining consistent naming conventions.
- It supports the hyperbolic polar form across all sectors of the plane.
- The matrix representation feature is based upon the robust foundation of [nalgebra::Matrix](https://docs.rs/nalgebra/latest/nalgebra/base/struct.Matrix.html).

## Usage

## Installation
`cargo add perplex_num` or add this to `Cargo.toml`:

```toml
[dependencies]
perplex_num = "0.1"
```

The `matrix` feature is enabled by default, which adds the `nalgebra` crate as a dependency. This can be disabled with:
```toml
[dependencies.perplex_num]
perplex_num = "0.1"
default-features = false
```
## Examples

The `examples` directory contains various practical demonstrations of how to use the `perplex_num` crate. These examples not only illustrate the usage of perplex numbers but also show how to produce visualizations as seen in the [Perplex Number Description](https://github.com/tomtuamnuq/perplex_num/blob/main/Perplex.md).

For instance, `examples/visualize_functions.rs` is executed by the following command:

```sh
cargo run --example visualize_functions
```

This will generate an image that depicts the behavior of functions like `sinh`, `cos`, `inv`, and `exp` when applied to perplex numbers.

### Creating a Perplex Number and Performing Operations

Here's a quick example of how to get started with creating a perplex number and performing basic operations:

```rust
use perplex_num::Perplex;

fn main() {
    // Create a Perplex number with real part 1.0 and hyperbolic part 0.5
    let z = Perplex::new(1.0, 0.5);
    // Calculate the hyperbolic sine of the perplex number
    let z_sinh = z.sinh();
    // Raise the perplex number or it's inverse to the power of 2
    let z_squared = z.powu(2);
    let z_inv_squared = z.powi(-2).expect("z is invertible");

    println!("The hyperbolic sine of {} is {:.5}", z, z_sinh);
    println!("{} raised to the power of 2 is {:.3}", z, z_squared);
    println!("{} raised to the power of -2 is {:.3}", z, z_inv_squared);
}
```

## Coverage

Test coverage report is generated with [cargo tarpaulin](https://github.com/xd009642/tarpaulin). Invoke it with:
```sh
cargo tarpaulin --verbose --all-targets --skip-clean --exclude-files "examples/*.rs" "benches/*.rs"
```


## Compatibility
The `perplex_num` crate is tested for rustc 1.76.

## Bibliography
- [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6)
- [Hyperbolic trigonometry in two-dimensional space-time geometry](https://doi.org/10.1393/ncb/i2003-10012-9)
- [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643)
- [Introduction to Hybrid Numbers](https://doi.org/10.1007/s00006-018-0833-3)
- [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586)
