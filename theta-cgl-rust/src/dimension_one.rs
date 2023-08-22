#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element
macro_rules! define_dim_one_theta_core{ () => {
    // Compute the Hadamard transform
    fn to_hadamard(X: Fq, Z: Fq) -> (Fq, Fq) {
        let X_new = X + Z;
        let Z_new = X - Z;

        (X_new, Z_new)
    }

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

        // Squared theta first squares the coords
        // then returns the hadamard transform. 
        // This gives the square of the dual coords
        pub fn squared_theta(self) -> (Fq, Fq) {
            let XX = self.X.square();
            let ZZ = self.Z.square();

            to_hadamard(XX, ZZ)

        }

        // Compute the two isogeny
        pub fn radical_two_isogeny(self, bit: u8) -> ThetaPoint {
            let (AA, BB) = self.squared_theta();
            let AABB = AA * BB; 
            let (mut AB, _) = AABB.sqrt();

            // TODO: make constant time
            let ctl = ((bit as u32) & 1).wrapping_neg();
            AB.set_condneg(ctl);

            let (X_new, Z_new) = to_hadamard(AA, AB);

            Self {
                X : X_new,
                Z : Z_new
            }
        }

        pub fn to_hash(self) -> Fq{
            self.Z / self.X
        }
    }

    pub fn cgl_hash(O0: ThetaPoint, msg: &[u8]) -> Fq{
        let mut r = O0;
        for bit in msg{
            r = r.radical_two_isogeny(*bit)
        }
        r.to_hash()
    }
} } // End of macro: define_dim_one_theta_core

pub(crate) use define_dim_one_theta_core;