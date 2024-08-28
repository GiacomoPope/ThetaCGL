#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element Fp^2
// X0, Z0, U0, V0, type Fq, coordinates of theta null point
macro_rules! define_dim_two_theta_core {
    () => {
        use crate::util::pad_msg;

        // Theta Point
        // The domain / codomain is described by a theta point
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaPointDim2 {
            pub X: Fq,
            pub Z: Fq,
            pub U: Fq,
            pub V: Fq,
        }

        impl ThetaPointDim2 {
            pub const fn new(X: &Fq, Z: &Fq, U: &Fq, V: &Fq) -> ThetaPointDim2 {
                Self {
                    X: *X,
                    Z: *Z,
                    U: *U,
                    V: *V,
                }
            }

            pub fn coords(self) -> (Fq, Fq, Fq, Fq) {
                (self.X, self.Z, self.U, self.V)
            }

            // Compute the Hadamard transform
            // Cost 8 additions
            fn to_hadamard(self, X: &Fq, Z: &Fq, U: &Fq, V: &Fq) -> (Fq, Fq, Fq, Fq) {
                let t1 = X + Z;
                let t2 = X - Z;
                let t3 = U + V;
                let t4 = U - V;

                (&t1 + &t3, &t2 + &t4, &t1 - &t3, &t2 - &t4)
            }

            // Squared theta first squares the coords
            // then returns the hadamard transform.
            // This gives the square of the dual coords
            pub fn squared_theta(self) -> (Fq, Fq, Fq, Fq) {
                let XX = self.X.square();
                let ZZ = self.Z.square();
                let UU = self.U.square();
                let VV = self.V.square();

                (XX, ZZ, UU, VV)
            }

            // Compute the two isogeny
            pub fn radical_two_isogeny(self, bits: Vec<u8>) -> ThetaPointDim2 {
                let (mut AA, mut BB, mut CC, mut DD) = self.squared_theta();

                (AA, BB, CC, DD) = self.to_hadamard(&AA, &BB, &CC, &DD);

                let AABB = &AA * &BB;
                let AACC = &AA * &CC;
                let AADD = &AA * &DD;

                let mut AB = AABB.sqrt().0;
                let mut AC = AACC.sqrt().0;
                let mut AD = AADD.sqrt().0;

                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                AB.set_condneg(ctl1);

                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                AC.set_condneg(ctl2);

                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                AD.set_condneg(ctl3);

                let (anew, bnew, cnew, dnew) = self.to_hadamard(&AA, &AB, &AC, &AD);

                ThetaPointDim2::new(&anew, &bnew, &cnew, &dnew)
            }

            // Compute a four-radical isogeny
            pub fn radical_four_isogeny(self, bits: Vec<u8>) -> ThetaPointDim2 {
                let (mut AA, mut BB, mut CC, mut DD) = self.squared_theta();
                (AA, BB, CC, DD) = self.to_hadamard(&AA, &BB, &CC, &DD);

                let AABB = &AA * &BB;
                let AACC = &AA * &CC;
                let BBDD = &BB * &DD;
                let CCDD = &CC * &DD;

                // Compute ABCD and flip the sign if bit[0] is set to 1
                let mut ABCD = (&AABB * &CCDD).sqrt().0;
                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                ABCD.set_condneg(ctl1);

                // Compute alpha_1 and alpha2
                let alpha_1_4 = (&ABCD.mul2() + &AABB + &CCDD).mul4();
                let alpha_2_4 = (&ABCD.mul2() + &AACC + &BBDD).mul4();
                let mut alpha_1 = alpha_1_4.fourth_root().0;
                let mut alpha_2 = alpha_2_4.fourth_root().0;

                // Multiply alpha_1 by a fourth root of unity
                // depending on bits 1 and 2
                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                alpha_1.set_cond(&(&alpha_1 * &Fq::ZETA), ctl3);
                alpha_1.set_condneg(ctl2);

                // Multiply alpha_2 by a fourth root of unity
                // depending on bits 3 and 4
                let ctl4 = ((bits[3] as u32) & 1).wrapping_neg();
                let ctl5 = ((bits[4] as u32) & 1).wrapping_neg();
                alpha_2.set_cond(&(&alpha_2 * &Fq::ZETA), ctl5);
                alpha_2.set_condneg(ctl4);

                let mut alpha_3_2 = (&CCDD + &ABCD).mul8();
                alpha_3_2 *= (&AACC + &ABCD) * &CCDD * &DD + &(&BBDD + &ABCD) * &CCDD * &CC;
                let mut alpha_3 = alpha_3_2.sqrt().0;

                // Change the sign of alpha_3 depending on bit 5
                let ctl6 = ((bits[5] as u32) & 1).wrapping_neg();
                alpha_3.set_condneg(ctl6);

                let lambda = &CCDD * &alpha_1 * &alpha_2;

                let (anew, bnew, cnew, dnew) = self.to_hadamard(
                    &(&self.X.mul2() * &lambda),
                    &(&alpha_1 * &lambda),
                    &(&alpha_2 * &lambda),
                    &alpha_3,
                );

                ThetaPointDim2::new(&anew, &bnew, &cnew, &dnew)
            }

            pub fn to_hash(self) -> (Fq, Fq, Fq) {
                let (X, Z, U, V) = (self.X, self.Z, self.U, self.V);
                let X_inv = X.invert();

                (&Z * &X_inv, &U * &X_inv, &V * &X_inv)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim2Rad2 {}

        impl CGLDim2Rad2 {
            const O0: ThetaPointDim2 = ThetaPointDim2::new(&X0, &Z0, &U0, &V0);

            pub fn new() -> Self {
                Self {}
            }
            pub fn bit_string(self, mut T: ThetaPointDim2, mut msg: Vec<u8>) -> ThetaPointDim2 {
                let chunk_len = 3;
                msg = pad_msg(msg, chunk_len);
                let iter = msg.chunks(chunk_len);
                for i in iter {
                    T = T.radical_two_isogeny(i.to_vec());
                }

                T
            }

            pub fn hash(self, msg: Vec<u8>) -> (Fq, Fq, Fq) {
                let T = self.bit_string(Self::O0, msg);
                T.to_hash()
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim2Rad4 {}

        impl CGLDim2Rad4 {
            const O0: ThetaPointDim2 = ThetaPointDim2::new(&X0, &Z0, &U0, &V0);

            pub fn new() -> Self {
                Self {}
            }

            pub fn bit_string(self, mut T: ThetaPointDim2, mut msg: Vec<u8>) -> ThetaPointDim2 {
                let chunk_len = 6;
                msg = pad_msg(msg, chunk_len);
                let iter = msg.chunks(chunk_len);
                for i in iter {
                    T = T.radical_four_isogeny(i.to_vec());
                }

                T
            }

            pub fn hash(&self, msg: Vec<u8>) -> (Fq, Fq, Fq) {
                let T = self.bit_string(Self::O0, msg);
                T.to_hash()
            }
        }
    };
} // End of macro: define_dim_two_theta_core

pub(crate) use define_dim_two_theta_core;
