All $\mathbb{F}_p$ and $\mathbb{F}_p^2$ arithmetic is created from the macro `fpcore.rs`. 
To generate the constants for the params, run `gen_fp.sage`. 

Code in `fpcore.rs` was written by Thomas Pornin for ongoing work in another project. I have
asked his permission to use this code here.

All arithmetic is implemented in constant time.

Test the arithmetic with
```
cargo test
```

Get benchmarks with

```
cargo bench
```
