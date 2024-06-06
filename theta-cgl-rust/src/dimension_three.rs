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
            pub const fn new(
                X: &Fq,
                Z: &Fq,
                U: &Fq,
                V: &Fq,
                G: &Fq,
                H: &Fq,
                I: &Fq,
                J: &Fq,
            ) -> ThetaPointDim3 {
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
                (
                    self.X, self.Z, self.U, self.V, self.G, self.H, self.I, self.J,
                )
            }

            // Compute the Hadamard transform
            fn to_hadamard(
                self,
                x0: Fq,
                x1: Fq,
                x2: Fq,
                x3: Fq,
                x4: Fq,
                x5: Fq,
                x6: Fq,
                x7: Fq,
            ) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
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

                let aa = a.square();
                let bb = b.square();
                let cc = c.square();
                let dd = d.square();
                let ee = e.square();
                let ff = f.square();
                let gg = g.square();
                let hh = h.square();

                let (mut AA, BB, CC, DD, EE, FF, GG, HH) =
                    self.to_hadamard(aa, bb, cc, dd, ee, ff, gg, hh);

                let mut lam = Fq::ZERO;

                let AA_is_non_zero = AA.iszero() ^ 0xFFFFFFFF;
                lam.set_cond(&AA, AA_is_non_zero);
                let mut non_zero_found = AA_is_non_zero;

                let mut all_non_zero = AA_is_non_zero;

                let mut set_lam = |x: Fq| {
                    let x_is_zero = x.iszero();
                    let x_is_non_zero = x_is_zero ^ 0xFFFFFFFF;
                    let non_zero_not_found = non_zero_found ^ 0xFFFFFFFF;
                    lam.set_cond(&x, non_zero_not_found & x_is_non_zero);
                    non_zero_found = non_zero_found | x_is_non_zero;
                    all_non_zero = all_non_zero & x_is_non_zero;
                };

                set_lam(BB);
                set_lam(CC);
                set_lam(DD);
                set_lam(EE);
                set_lam(FF);
                set_lam(GG);
                set_lam(HH);

                assert!(non_zero_found == 0xFFFFFFFF);

                let AAAA = lam * AA;
                let AABB = lam * BB;
                let AACC = lam * CC;
                let AADD = lam * DD;
                let AAEE = lam * EE;
                let AAFF = lam * FF;
                let AAGG = lam * GG;
                let AAHH = lam * HH;

                let mut AB = AABB.sqrt().0;
                let mut AC = AACC.sqrt().0;
                let mut AD = AADD.sqrt().0;
                let mut AE = AAEE.sqrt().0;
                let mut AF = AAFF.sqrt().0;
                let mut AG = AAGG.sqrt().0;

                // TODO
                // This seems to be unused?

                // let AB_is_zero = AB.iszero();
                // let mut is_some_zero = AB_is_zero;
                // is_some_zero = is_some_zero | AC.iszero();
                // is_some_zero = is_some_zero | AD.iszero();
                // is_some_zero = is_some_zero | AE.iszero();
                // is_some_zero = is_some_zero | AF.iszero();
                // is_some_zero = is_some_zero | AG.iszero();

                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                AB.set_condneg(ctl1);

                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                AC.set_condneg(ctl2);

                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                AD.set_condneg(ctl3);

                let ctl4 = ((bits[3] as u32) & 1).wrapping_neg();
                AE.set_condneg(ctl4);

                let ctl5 = ((bits[4] as u32) & 1).wrapping_neg();
                AF.set_condneg(ctl5);

                let ctl6 = ((bits[5] as u32) & 1).wrapping_neg();
                AG.set_condneg(ctl6);

                let AH;
                (AA, AB, AC, AD, AE, AF, AG, AH) = self.last_sqrt(
                    AAAA,
                    AABB,
                    AACC,
                    AADD,
                    AAEE,
                    AAFF,
                    AAGG,
                    AAHH,
                    AA,
                    AB,
                    AC,
                    AD,
                    AE,
                    AF,
                    AG,
                    lam,
                    all_non_zero,
                );

                let (anew, bnew, cnew, dnew, enew, fnew, gnew, hnew) =
                    self.to_hadamard(AA, AB, AC, AD, AE, AF, AG, AH);

                ThetaPointDim3::new(&anew, &bnew, &cnew, &dnew, &enew, &fnew, &gnew, &hnew)
            }

            fn last_sqrt(
                self,
                c0: Fq,
                c1: Fq,
                c2: Fq,
                c3: Fq,
                c4: Fq,
                c5: Fq,
                c6: Fq,
                c7: Fq,
                mut x0: Fq,
                mut x1: Fq,
                mut x2: Fq,
                mut x3: Fq,
                mut x4: Fq,
                mut x5: Fq,
                mut x6: Fq,
                lam: Fq,
                all_non_zero: u32,
            ) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let (a0, a1, a2, a3, a4, a5, a6, a7) =
                    self.to_hadamard(c0, c1, c2, c3, c4, c5, c6, c7);

                let R1 = a0 * a1 * a2 * a3;
                let R3 = a4 * a5 * a6 * a7;

                let c04 = c0 * c4;
                let c15 = c1 * c5;
                let c26 = c2 * c6;
                let c37 = c3 * c7;
                let c0246 = c04 * c26;
                let c1357 = c15 * c37;

                // TODO: implement mul_small in fp2_64_gen
                let mut tmp = (c04 - c15 + c26 - c37);
                tmp = tmp * tmp;
                let term = tmp.mul8().mul2() - (c0246 + c1357).mul8().mul8();

                tmp = (R1 + R3 - term);
                tmp = tmp * tmp;
                let num =
                    tmp + (c0246 * c1357).mul8().mul8().mul8().mul8().mul4() - (R1 * R3).mul4();
                let den = (R1 + R3 - term).mul8().mul8().mul4() * x0 * x1 * x2 * x3 * x4 * x5 * x6;

                let (y0, y1, y2, y3, y4, y5, y6, y7) = self.coords();
                let yy = y0 * y1 * y2 * y3 * y4 * y5 * y6 * y7;

                let lam_to_4 = lam * lam * lam * lam;
                let mut l = lam_to_4 * yy;
                l = l.mul8().mul8();
                let cnd1 = (c0246 * c1357).equals(&(l * l));

                let mut x7 = c7.sqrt().0;
                let den_is_zero = den.iszero();

                let tmp = -l / (x0 * x1 * x2 * x3 * x4 * x5 * x6);
                x7.set_cond(&tmp, den_is_zero & all_non_zero & cnd1);
                let p = x0 * x1 * x2 * x3 * x4 * x5 * x6;
                x7.set_cond(&(-l), den_is_zero & all_non_zero & cnd1);
                x0.set_cond(&(p * x0), den_is_zero & all_non_zero & cnd1);
                x1.set_cond(&(p * x1), den_is_zero & all_non_zero & cnd1);
                x2.set_cond(&(p * x2), den_is_zero & all_non_zero & cnd1);
                x3.set_cond(&(p * x3), den_is_zero & all_non_zero & cnd1);
                x4.set_cond(&(p * x4), den_is_zero & all_non_zero & cnd1);
                x5.set_cond(&(p * x5), den_is_zero & all_non_zero & cnd1);
                x6.set_cond(&(p * x6), den_is_zero & all_non_zero & cnd1);

                let den_is_not_zero = den_is_zero ^ 0xFFFFFFFF;
                x0.set_cond(&(den * x0), den_is_not_zero);
                x1.set_cond(&(den * x1), den_is_not_zero);
                x2.set_cond(&(den * x2), den_is_not_zero);
                x3.set_cond(&(den * x3), den_is_not_zero);
                x4.set_cond(&(den * x4), den_is_not_zero);
                x5.set_cond(&(den * x5), den_is_not_zero);
                x6.set_cond(&(den * x6), den_is_not_zero);
                x7.set_cond(&num, den_is_not_zero);

                (x0, x1, x2, x3, x4, x5, x6, x7)
            }

            pub fn to_hash(self) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let (X, Z, U, V, G, H, I, J) = (
                    self.X, self.Z, self.U, self.V, self.G, self.H, self.I, self.J,
                );
                let X_inv = X.invert();

                (
                    Z * X_inv,
                    U * X_inv,
                    V * X_inv,
                    G * X_inv,
                    H * X_inv,
                    I * X_inv,
                    J * X_inv,
                )
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
