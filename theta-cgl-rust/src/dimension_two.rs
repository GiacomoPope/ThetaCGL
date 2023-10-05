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
            fn to_hadamard(self, X: Fq, Z: Fq, U: Fq, V: Fq) -> (Fq, Fq, Fq, Fq) {
                (X + Z + U + V, X - Z + U - V, X + Z - U - V, X - Z - U + V)
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

                (AA, BB, CC, DD) = self.to_hadamard(AA, BB, CC, DD);

                let AABB = AA * BB;
                let AACC = AA * CC;
                let AADD = AA * DD;
                let mut AB = AABB.sqrt().0;
                let mut AC = AACC.sqrt().0;
                let mut AD = AADD.sqrt().0;

                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                AB.set_condneg(ctl1);

                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                AC.set_condneg(ctl2);

                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                AD.set_condneg(ctl3);

                let (anew, bnew, cnew, dnew) = self.to_hadamard(AA, AB, AC, AD);
                ThetaPointDim2::new(&anew, &bnew, &cnew, &dnew)
            }

            pub fn to_hash(self) -> (Fq, Fq, Fq) {
                let (X, Z, U, V) = (self.X, self.Z, self.U, self.V);
                let X_inv = X.invert();

                (Z * X_inv, U * X_inv, V * X_inv)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim2Rad2 {}

        impl CGLDim2Rad2 {
            const O0: ThetaPointDim2 = ThetaPointDim2::new(&X0, &Z0, &U0, &V0);

            pub fn new() -> Self {
                Self {}
            }
            pub fn bit_string(
                self,
                mut T: ThetaPointDim2,
                mut msg: Vec<u8>,
                chunk_len: usize,
            ) -> ThetaPointDim2 {
                msg = pad_msg(msg, chunk_len);
                let iter = msg.chunks(chunk_len);
                for i in iter {
                    T = T.radical_two_isogeny(i.to_vec());
                }

                T
            }

            pub fn hash(self, msg: Vec<u8>, chunk_len: usize) -> (Fq, Fq, Fq) {
                let T = self.bit_string(Self::O0, msg, chunk_len);
                T.to_hash()
            }
        }
    };
} // End of macro: define_dim_two_theta_core

pub(crate) use define_dim_two_theta_core;
