#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element
macro_rules! define_dim_two_theta_core{ () => {
    // Compute the Hadamard transform
    fn to_hadamard_dim2(X: Fq, Z: Fq, U: Fq, V: Fq) -> (Fq, Fq, Fq, Fq) {
        (X + Z + U + V,  X - Z + U - V, X + Z - U - V, X - Z - U + V)
    }

    // Theta Point
    // The domain / codomain is described by a theta point
    #[derive(Clone, Copy, Debug)]
    pub struct ThetaPointDim2 {
        pub X : Fq,
        pub Z : Fq,
        pub U : Fq,
        pub V : Fq
    }

    impl ThetaPointDim2 {    

        pub fn new(X: &Fq, Z: &Fq, U: &Fq, V: &Fq) -> ThetaPointDim2 {
            Self{X: *X, Z: *Z, U: *U, V: *V}
        }  

        pub fn coords(self) -> (Fq, Fq, Fq, Fq) {
            (self.X, self.Z, self.U, self.V)
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

            (AA, BB, CC, DD) = to_hadamard_dim2(AA, BB, CC, DD);

            let AABB = AA * BB;
            let AACC = AA * CC;
            let AADD = AA * DD;
            let mut AB = AABB.sqrt().0;
            let mut AC = AACC.sqrt().0;
            let mut AD = AADD.sqrt().0;
            
            if (bits[0] == 1) {
                AB = - AB
            }
            if (bits[1] == 1){
                AC = - AC
            }
            if (bits[2] == 1){
                AD = - AD
            }

            let (anew, bnew, cnew, dnew) = to_hadamard_dim2(AA, AB, AC, AD);
            ThetaPointDim2::new(&anew, &bnew, &cnew, &dnew)
        }

        pub fn to_hash(self) -> (Fq, Fq, Fq) {
            let (X, Z, U, V) = (self.X, self.Z, self.U, self.V);
            let X_inv = X.invert();

            (Z * X_inv, U * X_inv, V * X_inv)
        } 
    }

    pub fn bit_string(mut T: ThetaPointDim2, mut msg: Vec<u8>, chunk_len: usize) -> ThetaPointDim2 {
        let m = msg.len() % chunk_len;
        if m != 0 {
            for _ in 0..(chunk_len - m) {
                msg.push(0);
            }
        }
        let iter = msg.chunks(chunk_len);
        for i in iter {
            T = T.radical_two_isogeny(i.to_vec());
        }

        T
    }

    pub fn cgl_hash_dim2(O0: ThetaPointDim2, msg: Vec<u8>) -> (Fq, Fq, Fq) {
        let T = bit_string(O0, msg, 3);

        T.to_hash()
    }
} } // End of macro: define_dim_two_theta_core

pub(crate) use define_dim_two_theta_core;
