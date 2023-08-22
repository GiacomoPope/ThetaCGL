#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element
macro_rules! define_dim_one_theta_core{ () => {
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

        pub fn to_hash(self) -> Fq{
            self.Z / self.X
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CGL {
        pub O0 : ThetaPoint,
    }

    impl CGL {    

        pub fn new(O0: ThetaPoint) -> CGL {
            Self{O0: O0}
        }

        pub fn bit_string(self, mut T: ThetaPoint, msg: Vec<u8>) -> ThetaPoint {
            for bit in msg {
                T = T.radical_two_isogeny(bit)
            }

            T
        }

        pub fn hash(self, msg: Vec<u8>) -> Fq {
            let T = self.bit_string(self.O0, msg);

            T.to_hash()
        }
    }

} } // End of macro: define_dim_one_theta_core

pub(crate) use define_dim_one_theta_core;
