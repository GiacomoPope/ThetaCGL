#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element
// X0, Z0, U0, V0, type Fq, coordinates of theta null point
macro_rules! define_dim_one_theta_core {
    () => {
        use crate::util::pad_msg;

        // Theta Point
        // The domain / codomain is described by a theta point
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaPointDim1 {
            pub X: Fq,
            pub Z: Fq,
        }

        impl ThetaPointDim1 {
            pub const fn new(X: &Fq, Z: &Fq) -> ThetaPointDim1 {
                Self { X: *X, Z: *Z }
            }

            pub fn coords(self) -> (Fq, Fq) {
                (self.X, self.Z)
            }

            // Compute the Hadamard transform
            fn to_hadamard(self, X: &Fq, Z: &Fq) -> (Fq, Fq) {
                let X_new = X + Z;
                let Z_new = X - Z;

                (X_new, Z_new)
            }

            // Squared theta first squares the coords
            // then returns the hadamard transform.
            // This gives the square of the dual coords
            pub fn squared_theta(self) -> (Fq, Fq) {
                let XX = self.X.square();
                let ZZ = self.Z.square();

                self.to_hadamard(&XX, &ZZ)
            }

            // Compute the two isogeny
            pub fn radical_two_isogeny(self, bit: u8) -> ThetaPointDim1 {
                let (x0, x1) = self.squared_theta();
                let x01 = &x0 * &x1;
                let (mut y1, _) = x01.sqrt();

                let ctl = ((bit as u32) & 1).wrapping_neg();
                y1.set_condneg(ctl);

                let (b0, b1) = self.to_hadamard(&x0, &y1);

                Self { X: b0, Z: b1 }
            }

            pub fn radical_four_isogeny(self, bits: Vec<u8>) -> ThetaPointDim1 {
                let (x0, x1) = self.squared_theta();
                let x01 = &x0 * &x1;
                let mut factor = x01.fourth_root().0;

                // if the second bit is zero, we multiply by zeta
                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                factor.set_cond(&(&factor * &Fq::ZETA), ctl2);

                // if the first bit is zero, we negate the result
                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                factor.set_condneg(ctl1);

                let b0 = self.X + factor;
                let b1 = self.X - factor;
                let (b0, b1) = self.to_hadamard(&b0, &b1);

                Self { X: b0, Z: b1 }
            }

            pub fn radical_eight_isogeny(
                self,
                bits: Vec<u8>,
                torsion: ThetaPointDim1,
                zeta_8: &Fq,
                sqrt_two: &Fq,
            ) -> (ThetaPointDim1, ThetaPointDim1) {
                let (a0, a1) = self.coords();
                let (u0, u1) = torsion.coords();

                let a00 = a0.square();
                let a01 = &a0 * &a1;
                let a11 = a1.square();

                let u00 = u0.square();
                let u01 = &u0 * &u1;

                let u0_8 = u00.square().square();
                let u1_8 = u1.square().square().square();

                let mut factor = (&u0_8 - &u1_8).eighth_root().0;

                // if the first bit is zero, we negate the result
                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                factor.set_condneg(ctl1);

                // if the second bit is zero, we multiply by 4th root of unity
                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                factor.set_cond(&(&factor * &Fq::ZETA), ctl2);

                // if the third bit is zero, we multiply by 8th root of unity
                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                factor.set_cond(&(&factor * zeta_8), ctl3);

                let factor_2 = factor.square();
                let factor_4 = factor_2.square();

                // Compute the codomain
                let (b0, b1) = self.to_hadamard(&u00, &factor_2);
                let B = Self { X: b0, Z: b1 };

                // Reconstruct the torsion for the next step
                let t = &a01 * &u01.mul2();
                let v0 = &t * &b1;
                let v1 =
                    &a00.mul2() * &u01.square() + &factor_4 * &a11 - sqrt_two * &factor * &t * &u0;
                let T = Self { X: v0, Z: v1 };

                (B, T)
            }

            pub fn to_hash(self) -> Fq {
                self.Z / self.X
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim1Rad2 {
            block_size: usize,
        }

        impl CGLDim1Rad2 {
            const O0: ThetaPointDim1 = ThetaPointDim1::new(&X0, &Z0);

            pub fn new(block_size: usize) -> CGLDim1Rad2 {
                Self { block_size }
            }

            pub fn bit_string(&self, mut T: ThetaPointDim1, msg: Vec<u8>) -> ThetaPointDim1 {
                for bit in msg {
                    T = T.radical_two_isogeny(bit)
                }

                T
            }

            pub fn hash(&self, msg: Vec<u8>) -> Fq {
                let padded_msg = pad_msg(msg, self.block_size);
                let T = self.bit_string(Self::O0, padded_msg);
                T.to_hash()
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim1Rad4 {
            chunk_len: usize,
            block_size: usize,
        }

        impl CGLDim1Rad4 {
            const O0: ThetaPointDim1 = ThetaPointDim1::new(&X0, &Z0);

            pub fn new(block_size: usize) -> CGLDim1Rad4 {
                let chunk_len = 2;
                assert!(block_size % chunk_len == 0);
                Self {
                    chunk_len,
                    block_size,
                }
            }

            pub fn bit_string(self, mut T: ThetaPointDim1, msg: Vec<u8>) -> ThetaPointDim1 {
                let iter = msg.chunks(self.chunk_len);
                for i in iter {
                    T = T.radical_four_isogeny(i.to_vec());
                }

                T
            }

            pub fn hash(&self, msg: Vec<u8>) -> Fq {
                let padded_msg = pad_msg(msg, self.block_size);
                let T = self.bit_string(Self::O0, padded_msg);
                T.to_hash()
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim1Rad8 {
            chunk_len: usize,
            block_size: usize,
        }

        impl CGLDim1Rad8 {
            // Null point
            const O0: ThetaPointDim1 = ThetaPointDim1::new(&X0, &Z0);

            // Torsion point
            const TORSION: ThetaPointDim1 = ThetaPointDim1::new(&PX0, &PZ0);

            // Sqrt of Two
            const SQRT_TWO: Fq = fp2_sqrt_2;

            // eighth root of unity
            const ZETA_8: Fq = fp2_zeta_8;

            pub fn new(block_size: usize) -> CGLDim1Rad8 {
                let chunk_len = 3;
                assert!(block_size % chunk_len == 0);
                Self {
                    chunk_len,
                    block_size,
                }
            }

            pub fn bit_string(self, mut T: ThetaPointDim1, msg: Vec<u8>) -> ThetaPointDim1 {
                let iter = msg.chunks(self.chunk_len);

                // Set the torsion point
                let mut torsion = Self::TORSION;
                for i in iter {
                    (T, torsion) = T.radical_eight_isogeny(
                        i.to_vec(),
                        torsion,
                        &Self::ZETA_8,
                        &Self::SQRT_TWO,
                    );
                }

                T
            }

            pub fn hash(&self, msg: Vec<u8>) -> Fq {
                let padded_msg = pad_msg(msg, self.block_size);
                let T = self.bit_string(Self::O0, padded_msg);
                T.to_hash()
            }
        }
    };
} // End of macro: define_dim_one_theta_core

pub(crate) use define_dim_one_theta_core;
