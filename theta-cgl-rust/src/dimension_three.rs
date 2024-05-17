#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element Fp^2
// X0, Z0, U0, V0, G0, H0, I0, J0 type Fq, coordinates of theta null point
macro_rules! define_dim_three_theta_core {
    () => {
        use crate::util::pad_msg;

        // Theta Point
        // The domain / codomain is described by a theta point
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaPointDim3 {
            pub X: Fq,
            pub Z: Fq,
            pub U: Fq,
            pub V: Fq,
            pub G: Fq,
            pub H: Fq,
            pub I: Fq,
            pub J: Fq,
        }

        impl ThetaPointDim3 {
            pub const fn new(X: &Fq, Z: &Fq, U: &Fq, V: &Fq, G: &Fq, H: &Fq, I: &Fq, J: &Fq) -> ThetaPointDim3 {
                Self {
                    X: *X,
                    Z: *Z,
                    U: *U,
                    V: *V,
                    G: *G,
                    H: *H,
                    I: *I,
                    J: *J,
                }
            }

            pub fn coords(self) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                (self.X, self.Z, self.U, self.V, self.G, self.H, self.I, self.J)
            }

            // Compute the Hadamard transform
            fn to_hadamard(self, x0: Fq, x1: Fq, x2: Fq, x3: Fq, x4: Fq, x5: Fq, x6: Fq, x7: Fq) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let y0 = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;
                let y1 = x0 - x1 + x2 - x3 + x4 - x5 + x6 - x7;
                let y2 = x0 + x1 - x2 - x3 + x4 + x5 - x6 - x7;
                let y3 = x0 - x1 - x2 + x3 + x4 - x5 - x6 + x7;
                let y4 = x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7;
                let y5 = x0 - x1 + x2 - x3 - x4 + x5 - x6 + x7;
                let y6 = x0 + x1 - x2 - x3 - x4 - x5 + x6 + x7;
                let y7 = x0 - x1 - x2 + x3 - x4 + x5 + x6 - x7;

                (y0, y1, y2, y3, y4, y5, y6, y7)
            }

            // Compute the two isogeny
            pub fn radical_two_isogeny(self, bits: Vec<u8>) -> ThetaPointDim3 {
                let (a, b, c, d, e, f, g, h) = self.coords();

                let aa = a*a;
                let bb = b*b;
                let cc = c*c;
                let dd = d*d;
                let ee = e*e;
                let ff = f*f;
                let gg = g*g;
                let hh = h*h;
                let (AA, BB, CC, DD, EE, FF, GG, HH) = self.to_hadamard(aa, bb, cc, dd, ee, ff, gg, hh);

                let AABB = AA * BB;
                let AACC = AA * CC;
                let AADD = AA * DD;
                let AAEE = AA * EE;
                let AAFF = AA * FF;
                let AAGG = AA * GG;
                let AAHH = AA * HH;
                let mut AB = AABB.sqrt().0;
                let mut AC = AACC.sqrt().0;
                let mut AD = AADD.sqrt().0;
                let mut AE = AAEE.sqrt().0;
                let mut AF = AAFF.sqrt().0;
                let mut AG = AAGG.sqrt().0;

                let ctl1 = ((bits[0] as u64) & 1).wrapping_neg();
                AB.set_condneg(ctl1);

                let ctl2 = ((bits[1] as u64) & 1).wrapping_neg();
                AC.set_condneg(ctl2);

                let ctl3 = ((bits[2] as u64) & 1).wrapping_neg();
                AD.set_condneg(ctl3);

                let ctl4 = ((bits[3] as u64) & 1).wrapping_neg();
                AE.set_condneg(ctl4);

                let ctl5 = ((bits[4] as u64) & 1).wrapping_neg();
                AF.set_condneg(ctl5);

                let ctl6 = ((bits[5] as u64) & 1).wrapping_neg();
                AG.set_condneg(ctl6);

                // TODO:
                // AH =  ThetaCGLDim3.last_sqrt(self,AA*AA,AABB,AACC,AADD,AAEE,AAFF,AAGG,AAHH,AA,AB,AC,AD,AE,AF,AG)
                let AH = AAHH.sqrt().0;

                let (anew, bnew, cnew, dnew, enew, fnew, gnew, hnew) = self.to_hadamard(AA, AB, AC, AD, AE, AF, AG, AH);

                ThetaPointDim3::new(&anew, &bnew, &cnew, &dnew, &enew, &fnew, &gnew, &hnew)
            }

            pub fn to_hash(self) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let (X, Z, U, V, G, H, I, J) = (self.X, self.Z, self.U, self.V, self.G, self.H, self.I, self.J);
                let X_inv = X.invert();

                (Z * X_inv, U * X_inv, V * X_inv, G * X_inv, H * X_inv, I * X_inv, J * X_inv)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim3Rad2 {}

        impl CGLDim3Rad2 {
            const O0: ThetaPointDim3 = ThetaPointDim3::new(&X0, &Z0, &U0, &V0, &G0, &H0, &I0, &J0);

            pub fn new() -> Self {
                Self {}
            }
            pub fn bit_string(self, mut T: ThetaPointDim3, mut msg: Vec<u8>) -> ThetaPointDim3 {
                let chunk_len = 6;
                msg = pad_msg(msg, chunk_len);
                let iter = msg.chunks(chunk_len);
                for i in iter {
                    T = T.radical_two_isogeny(i.to_vec());
                }

                T
            }

            pub fn hash(self, msg: Vec<u8>) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let T = self.bit_string(Self::O0, msg);
                T.to_hash()
            }
        }
    };
} // End of macro: define_dim_three_theta_core

pub(crate) use define_dim_three_theta_core;
