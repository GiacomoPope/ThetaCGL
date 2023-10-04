#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element
macro_rules! define_dim_one_theta_core{ () => {
    use crate::util::pad_msg;

    // Theta Point
    // The domain / codomain is described by a theta point
    #[derive(Clone, Copy, Debug)]
    pub struct ThetaPointDim1 {
        pub X : Fq,
        pub Z : Fq
    }

    impl ThetaPointDim1 {    

        pub fn new(X: &Fq, Z: &Fq) -> ThetaPointDim1 {
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
        pub fn radical_two_isogeny(self, bit: u8) -> ThetaPointDim1 {
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

        pub fn radical_four_isogeny(self, bits: Vec<u8>, zeta: Fq) -> ThetaPointDim1 {
            let (AA, BB) = self.squared_theta();
            let AABB = AA * BB; 

            let mut factor = AABB.fourth_root().0;
           
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
    pub struct CGLDim1Rad2 {
        pub O0 : ThetaPointDim1,
    }

    impl CGLDim1Rad2 {    

        pub fn new(O0: ThetaPointDim1) -> CGLDim1Rad2 {
            Self{
                O0: O0,
            }
        }

        pub fn bit_string(&self, mut T: ThetaPointDim1, msg: Vec<u8>) -> ThetaPointDim1 {
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
    pub struct CGLDim1Rad4 {
        pub O0 : ThetaPointDim1,
        pub zeta: Fq,
    }

    impl CGLDim1Rad4 {

        pub fn new(O0: ThetaPointDim1, zeta: Fq) -> CGLDim1Rad4 {
            Self{
                O0: O0,
                zeta: zeta,
            }
        }

        pub fn bit_string(self, mut T: ThetaPointDim1, mut msg: Vec<u8>) -> ThetaPointDim1 {
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
