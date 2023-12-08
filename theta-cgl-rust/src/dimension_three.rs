#![allow(non_snake_case)]

// Macro expectations:
// Fq      type of field element Fp^2
// X0, Z0, U0, V0, S0, T0 type Fq, coordinates of theta null point
macro_rules! define_dim_three_theta_core {
    () => {
        use crate::util::pad_msg;

        // Theta Point
        // The domain / codomain is described by a theta point
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaPointDim3 {
            pub A: Fq,
            pub B: Fq,
            pub C: Fq,
            pub D: Fq,
            pub E: Fq,
            pub F: Fq,
            pub G: Fq,
            pub H: Fq,
        }

        impl ThetaPointDim3 {
            pub const fn new(A: &Fq, B: &Fq, C: &Fq, D: &Fq, E: &Fq, F: &Fq, G: &Fq, H: &Fq) -> ThetaPointDim3 {
                Self {
                    A: *A,
                    B: *B,
                    C: *C,
                    D: *D,
                    E: *E,
                    F: *F,
                    G: *G,
                    H: *H,
                }
            }

            pub fn coords(self) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                (self.A, self.B, self.C, self.D, self.E, self.F, self.G, self.H)
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
            pub fn radical_two_isogeny(self, bits: Vec<u8>) -> ThetaPointDim2 {
                let (a, b, c, d, e, f, g, h) = self.coord();

                let aa = a*a
                bb = b*b
                cc = c*c
                dd = d*d
                ee = e*e
                ff = f*f
                gg = g*g
                hh = h*h
                AA, BB, CC, DD, EE, FF, GG, HH = ThetaCGLDim3.hadamard(aa,bb,cc,dd,ee,ff,gg,hh)

                AABB = AA * BB
                AACC = AA * CC
                AADD = AA * DD
                AAEE = AA * EE
                AAFF = AA * FF
                AAGG = AA * GG
                AAHH = AA * HH
                AB = self.sqrt(AABB)
                AC = self.sqrt(AACC)
                AD = self.sqrt(AADD)
                AE = self.sqrt(AAEE)
                AF = self.sqrt(AAFF)
                AG = self.sqrt(AAGG)        
                
                if bits[0] == 1:
                    AB = - AB
                if bits[1] == 1:
                    AC = - AC
                if bits[2] == 1:
                    AD = - AD
                if bits[3] == 1:
                    AE = - AE
                if bits[4] == 1:
                    AF = - AF
                if bits[5] == 1:
                    AG = - AG

                AH =  ThetaCGLDim3.last_sqrt(self,AA*AA,AABB,AACC,AADD,AAEE,AAFF,AAGG,AAHH,AA,AB,AC,AD,AE,AF,AG)

                anew, bnew, cnew, dnew, enew, fnew, gnew, hnew = ThetaCGLDim3.hadamard(AA, AB, AC, AD, AE, AF, AG, AH)
                O1 = ThetaNullPointDim3(anew, bnew, cnew, dnew, enew, fnew, gnew, hnew)
                return O1

                /*
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
                */
            }

            pub fn to_hash(self) -> (Fq, Fq, Fq) {
                let (X, Z, U, V) = (self.X, self.Z, self.U, self.V);
                let X_inv = X.invert();

                (Z * X_inv, U * X_inv, V * X_inv)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CGLDim3Rad2 {}

        impl CGLDim3Rad2 {
            const O0: ThetaPointDim3 = ThetaPointDim3::new(&X0, &Z0, &U0, &V0, &S0, &T0);

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

            pub fn hash(self, msg: Vec<u8>) -> (Fq, Fq, Fq) {
                let T = self.bit_string(Self::O0, msg);
                T.to_hash()
            }
        }
    };
} // End of macro: define_dim_three_theta_core

pub(crate) use define_dim_three_theta_core;
