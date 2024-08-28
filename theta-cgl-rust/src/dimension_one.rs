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
                let (AA, BB) = self.squared_theta();
                let AABB = &AA * &BB;
                let (mut AB, _) = AABB.sqrt();

                let ctl = ((bit as u32) & 1).wrapping_neg();
                AB.set_condneg(ctl);

                let (X_new, Z_new) = self.to_hadamard(&AA, &AB);

                Self { X: X_new, Z: Z_new }
            }

            pub fn radical_four_isogeny(self, bits: Vec<u8>) -> ThetaPointDim1 {
                let (AA, BB) = self.squared_theta();
                let AABB = &AA * &BB;
                let (mut factor, check) = AABB.fourth_root();
                if check == 0 {
                    panic!("Something has gone wrong with the isogeny chain.")
                }

                // if the second bit is zero, we multiply by zeta
                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                factor.set_cond(&(&factor * &Fq::ZETA), ctl2);

                // if the first bit is zero, we negate the result
                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                factor.set_condneg(ctl1);

                let X_new = self.X + factor;
                let Z_new = self.X - factor;
                let (X_new, Z_new) = self.to_hadamard(&X_new, &Z_new);

                Self { X: X_new, Z: Z_new }
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
    };
} // End of macro: define_dim_one_theta_core

pub(crate) use define_dim_one_theta_core;
