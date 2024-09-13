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
                let (a, b) = self.coords();
                let (r, s) = torsion.coords();

                let rr = r.square();
                let r8 = rr.square().square();
                let s8 = s.square().square().square();

                let (mut factor, check) = (&r8 - &s8).eighth_root();
                if (check == 0) {
                    panic!("Something went wrong!");
                }

                // if the first bit is zero, we negate the result
                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                factor.set_condneg(ctl1);

                // if the second bit is zero, we multiply by 4th root of unity
                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                factor.set_cond(&(&factor * &Fq::ZETA), ctl2);

                // if the third bit is zero, we multiply by 8th root of unity
                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                factor.set_cond(&(&factor * zeta_8), ctl3);

                // Precompute some constants
                let rs = &r * &s;
                let abrs_2 = &a * &b * &rs.mul2();
                let factor_2 = factor.square();
                let factor_4 = factor_2.square();

                // Compute the codomain
                let b0 = &rr + &factor_2;
                let b1 = &rr - &factor_2;
                let B = Self { X: b0, Z: b1 };

                // Reconstruct the torsion for the next step
                let r4 = &abrs_2 * &b1;
                let s4 = &a.square().mul2() * &rs.square() + &factor_4 * &b.square()
                    - sqrt_two * &factor * &abrs_2 * &r;
                let T = Self { X: r4, Z: s4 };

                (B, T)
            }

            pub fn to_hash(self) -> Fq {
                self.Z / self.X
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim1Rad2 {}

        impl CGLDim1Rad2 {
            const O0: ThetaPointDim1 = ThetaPointDim1::new(&X0, &Z0);

            pub fn new() -> CGLDim1Rad2 {
                Self {}
            }

            pub fn bit_string(&self, mut T: ThetaPointDim1, msg: Vec<u8>) -> ThetaPointDim1 {
                for bit in msg {
                    T = T.radical_two_isogeny(bit)
                }

                T
            }

            pub fn hash(&self, msg: Vec<u8>) -> Fq {
                let T = self.bit_string(Self::O0, msg);

                T.to_hash()
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim1Rad4 {}

        impl CGLDim1Rad4 {
            const O0: ThetaPointDim1 = ThetaPointDim1::new(&X0, &Z0);

            pub fn new() -> CGLDim1Rad4 {
                Self {}
            }

            pub fn bit_string(self, mut T: ThetaPointDim1, mut msg: Vec<u8>) -> ThetaPointDim1 {
                let chunk_len = 2;
                msg = pad_msg(msg, chunk_len);
                let iter = msg.chunks(chunk_len);
                for i in iter {
                    T = T.radical_four_isogeny(i.to_vec());
                }

                T
            }

            pub fn hash(&self, msg: Vec<u8>) -> Fq {
                let T = self.bit_string(Self::O0, msg);

                T.to_hash()
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim1Rad8 {}

        impl CGLDim1Rad8 {
            // Null point
            const O0: ThetaPointDim1 = ThetaPointDim1::new(&X0, &Z0);

            // Torsion point
            const TORSION: ThetaPointDim1 = ThetaPointDim1::new(&PX0, &PZ0);

            // Sqrt of Two
            const SQRT_TWO: Fq = fp2_sqrt_2;

            // eighth root of unity
            const ZETA_8: Fq = fp2_zeta_8;

            pub fn new() -> CGLDim1Rad8 {
                Self {}
            }

            pub fn bit_string(self, mut T: ThetaPointDim1, mut msg: Vec<u8>) -> ThetaPointDim1 {
                let chunk_len = 3;
                msg = pad_msg(msg, chunk_len);
                let iter = msg.chunks(chunk_len);

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
                let T = self.bit_string(Self::O0, msg);

                T.to_hash()
            }
        }
    };
} // End of macro: define_dim_one_theta_core

pub(crate) use define_dim_one_theta_core;
