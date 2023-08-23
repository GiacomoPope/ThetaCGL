#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element
macro_rules! define_dim_one_theta_core{ () => {
    use crate::util::pad_msg;

    // Theta Point
    // The domain / codomain is described by a theta point
    #[derive(Clone, Copy, Debug)]
    pub struct ThetaPoint {
        pub X : Fq,
        pub Z : Fq
    }

    impl ThetaPoint {    

        pub fn new(X: &Fq, Z: &Fq) -> ThetaPoint {
            Self{X: *X, Z: *Z}
        }  

        pub fn coords(self) -> (Fq, Fq) {
            (self.X, self.Z)
        }

        // Compute the Hadamard transform
        fn to_hadamard(self, X: Fq, Z: Fq) -> (Fq, Fq) {
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

            self.to_hadamard(XX, ZZ)

        }

        // Compute the two isogeny
        pub fn radical_two_isogeny(self, bit: u8) -> ThetaPoint {
            let (AA, BB) = self.squared_theta();
            let AABB = AA * BB; 
            let (mut AB, _) = AABB.sqrt();

            let ctl = ((bit as u32) & 1).wrapping_neg();
            AB.set_condneg(ctl);

            let (X_new, Z_new) = self.to_hadamard(AA, AB);

            Self {
                X : X_new,
                Z : Z_new
            }
        }

        pub fn radical_four_isogeny(self, bits: Vec<u8>, zeta: Fq) -> ThetaPoint {
            let (AA, BB) = self.squared_theta();
            let AABB = AA * BB; 
            let AB = AABB.sqrt().0;
            let mut factor = AB.sqrt().0;

            let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
            factor.set_condneg(ctl1);

            // TODO: constant time
            if (bits[1] == 1) {
                factor = zeta * factor;
            }

            let X_new = self.X + factor;
            let Z_new = self.X - factor;

            // TODO: currently no hadamard call - as in Python code
            // anew, bnew = ThetaCGL.hadamard(anew, bnew) # I think we need an hadamard?
            // let (X_new, Z_new) = self.to_hadamard(X_new, Z_new);

            Self {
                X : X_new,
                Z : Z_new
            }
        }

        pub fn to_hash(self) -> Fq {
            self.Z / self.X
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CGL_1_2 {
        pub O0 : ThetaPoint,
    }

    impl CGL_1_2 {    

        pub fn new(O0: ThetaPoint) -> CGL_1_2 {
            Self{
                O0: O0,
            }
        }

        pub fn bit_string(&self, mut T: ThetaPoint, msg: Vec<u8>) -> ThetaPoint {
            for bit in msg {
                T = T.radical_two_isogeny(bit)
            }

            T
        } 

        pub fn hash(&self, msg: Vec<u8>) -> Fq {
            let T = self.bit_string(self.O0, msg);

            T.to_hash()
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CGL_1_4 {
        pub O0 : ThetaPoint,
        pub zeta: Fq,
    }

    impl CGL_1_4 {

        pub fn new(O0: ThetaPoint, zeta: Fq) -> CGL_1_4 {
            Self{
                O0: O0,
                zeta: zeta,
            }
        }

        pub fn bit_string(self, mut T: ThetaPoint, mut msg: Vec<u8>) -> ThetaPoint {
            let chunk_len = 2;
            msg = pad_msg(msg, chunk_len); 
            let iter = msg.chunks(chunk_len);
            for i in iter {
                T = T.radical_four_isogeny(i.to_vec(), self.zeta);
            }

            T
        }

        pub fn hash(&self, msg: Vec<u8>) -> Fq {
            let T = self.bit_string(self.O0, msg);

            T.to_hash()
        }
    }

} } // End of macro: define_dim_one_theta_core

pub(crate) use define_dim_one_theta_core;
