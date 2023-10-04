macro_rules! define_fp_bench {
    () => {
        mod util;
        use util::core_cycles;

        fn mkfp() -> Fp {
            let mut buf = [0u8; (Fp::ENCODED_LENGTH + 7) & !7usize];
            for i in 0..(buf.len() >> 3) {
                buf[(i << 3)..((i + 1) << 3)].copy_from_slice(&core_cycles().to_be_bytes());
            }
            Fp::decode_reduce(&buf)
        }

        fn bench_fp_add() {
            let mut x = mkfp();
            let mut y = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x += &y;
                    y += &x;
                    x += &y;
                    y += &x;
                    x += &y;
                    y += &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) add:            {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0]
            );
        }

        fn bench_fp_sub() {
            let mut x = mkfp();
            let mut y = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x -= &y;
                    y -= &x;
                    x -= &y;
                    y -= &x;
                    x -= &y;
                    y -= &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) sub:            {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0]
            );
        }

        fn bench_fp_mul_small() {
            let mut x = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let k = core_cycles() as i32;
                let begin = core_cycles();
                for _ in 0..1000 {
                    x.set_mul_small(k);
                    x.set_mul_small(k);
                    x.set_mul_small(k);
                    x.set_mul_small(k);
                    x.set_mul_small(k);
                    x.set_mul_small(k);
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) mul_small:      {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0]
            );
        }

        fn bench_fp_mul() {
            let mut x = mkfp();
            let mut y = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x *= &y;
                    y *= &x;
                    x *= &y;
                    y *= &x;
                    x *= &y;
                    y *= &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) mul:            {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0]
            );
        }

        fn bench_fp_square() {
            let mut x = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x.set_square();
                    x.set_square();
                    x.set_square();
                    x.set_square();
                    x.set_square();
                    x.set_square();
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) square:         {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0]
            );
        }

        fn bench_fp_div() {
            let mut x = mkfp();
            let mut y = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..100 {
                    x /= &y;
                    y /= &x;
                    x /= &y;
                    y /= &x;
                    x /= &y;
                    y /= &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) div:            {:13.2}  ({})",
                (tt[4] as f64) / 600.0,
                x.encode()[0]
            );
        }

        fn bench_fp_legendre() {
            let mut x = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..600 {
                    let ls = x.legendre();
                    x.set_cond(&x.mul2(), (ls >> 1) as u32);
                    x.set_cond(&x.mul4(), (-ls >> 1) as u32);
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) legendre:       {:13.2}  ({})",
                (tt[4] as f64) / 600.0,
                x.encode()[0]
            );
        }

        fn bench_fp_sqrt() {
            let mut x = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..20 {
                    let (mut x2, r) = x.sqrt();
                    x2.set_cond(&Fp::ONE, !r);
                    x += &x2;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) sqrt:           {:13.2}  ({})",
                (tt[4] as f64) / 20.0,
                x.encode()[0]
            );
        }

        fn bench_fp_fourth_root() {
            let mut x = mkfp();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..20 {
                    let (mut x2, r) = x.fourth_root();
                    x2.set_cond(&Fp::ONE, !r);
                    x += &x2;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p) fourth root:    {:13.2}  ({})",
                (tt[4] as f64) / 20.0,
                x.encode()[0]
            );
        }

        fn mkfp2() -> Fp2 {
            let re = mkfp();
            let im = mkfp();
            Fp2::new(&re, &im)
        }

        fn bench_fp2_add() {
            let mut x = mkfp2();
            let mut y = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x += &y;
                    y += &x;
                    x += &y;
                    y += &x;
                    x += &y;
                    y += &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) add:          {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_sub() {
            let mut x = mkfp2();
            let mut y = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x -= &y;
                    y -= &x;
                    x -= &y;
                    y -= &x;
                    x -= &y;
                    y -= &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) sub:          {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_mul() {
            let mut x = mkfp2();
            let mut y = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..1000 {
                    x *= &y;
                    y *= &x;
                    x *= &y;
                    y *= &x;
                    x *= &y;
                    y *= &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) mul:          {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_square() {
            let mut x = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..6000 {
                    x.set_square();
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) square:       {:13.2}  ({})",
                (tt[4] as f64) / 6000.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_div() {
            let mut x = mkfp2();
            let mut y = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..100 {
                    x /= &y;
                    y /= &x;
                    x /= &y;
                    y /= &x;
                    x /= &y;
                    y /= &x;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) div:          {:13.2}  ({})",
                (tt[4] as f64) / 600.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_legendre() {
            let mut x = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..600 {
                    let ls = x.legendre();
                    x.set_cond(&x.mul2(), (ls >> 1) as u32);
                    x.set_cond(&x.mul4(), (-ls >> 1) as u32);
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) legendre:     {:13.2}  ({})",
                (tt[4] as f64) / 600.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_sqrt() {
            let mut x = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..10 {
                    let (mut x2, r) = x.sqrt();
                    x2.set_cond(&Fp2::ONE, !r);
                    x += &x2;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) sqrt:         {:13.2}  ({})",
                (tt[4] as f64) / 10.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn bench_fp2_fourth_root() {
            let mut x = mkfp2();
            let mut tt = [0; 10];
            for i in 0..10 {
                let begin = core_cycles();
                for _ in 0..10 {
                    let (mut x2, r) = x.fourth_root();
                    x2.set_cond(&Fp2::ONE, !r);
                    x += &x2;
                }
                let end = core_cycles();
                tt[i] = end.wrapping_sub(begin);
            }
            tt.sort();
            println!(
                "GF(p^2) fourth_root   {:13.2}  ({})",
                (tt[4] as f64) / 10.0,
                x.encode()[0] ^ x.encode()[Fp::ENCODED_LENGTH]
            );
        }

        fn main() {
            println!("### p = {:?}", FP_NAME);
            bench_fp_add();
            bench_fp_sub();
            bench_fp_mul_small();
            bench_fp_mul();
            bench_fp_square();
            bench_fp_div();
            bench_fp_legendre();
            bench_fp_sqrt();
            bench_fp_fourth_root();
            bench_fp2_add();
            bench_fp2_sub();
            bench_fp2_mul();
            bench_fp2_square();
            bench_fp2_div();
            bench_fp2_legendre();
            bench_fp2_sqrt();
            bench_fp2_fourth_root();
        }
    };
} // End of macro: define_fp_bench

pub(crate) use define_fp_bench;
