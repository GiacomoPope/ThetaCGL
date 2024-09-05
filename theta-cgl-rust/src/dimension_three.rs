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
            // 24 additions
            #[inline(always)]
            fn to_hadamard(
                self,
                x0: &Fq,
                x1: &Fq,
                x2: &Fq,
                x3: &Fq,
                x4: &Fq,
                x5: &Fq,
                x6: &Fq,
                x7: &Fq,
            ) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let t0 = x0 + x1;
                let t1 = x0 - x1;
                let t2 = x2 + x3;
                let t3 = x2 - x3;
                let t4 = x4 + x5;
                let t5 = x4 - x5;
                let t6 = x6 + x7;
                let t7 = x6 - x7;

                let s0 = &t0 + &t2;
                let s1 = &t0 - &t2;
                let s2 = &t1 + &t3;
                let s3 = &t1 - &t3;
                let s4 = &t4 + &t6;
                let s5 = &t4 - &t6;
                let s6 = &t5 + &t7;
                let s7 = &t5 - &t7;

                let y0 = &s0 + &s4;
                let y1 = &s2 + &s6;
                let y2 = &s1 + &s5;
                let y3 = &s3 + &s7;
                let y4 = &s0 - &s4;
                let y5 = &s2 - &s6;
                let y6 = &s1 - &s5;
                let y7 = &s3 - &s7;

                (y0, y1, y2, y3, y4, y5, y6, y7)
            }

            // Compute the two isogeny
            pub fn radical_two_isogeny(self, bits: Vec<u8>) -> ThetaPointDim3 {
                let (a0, a1, a2, a3, a4, a5, a6, a7) = self.coords();

                let mut x0 = a0.square();
                let mut x1 = a1.square();
                let mut x2 = a2.square();
                let mut x3 = a3.square();
                let mut x4 = a4.square();
                let mut x5 = a5.square();
                let mut x6 = a6.square();
                let mut x7 = a7.square();
                (x0, x1, x2, x3, x4, x5, x6, x7) =
                    self.to_hadamard(&x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7);

                // We must consider an edge-case where exactly one of AA, ..., HH is
                // zero. To simplify the implementation we force the zero coordinate
                // to be HH temporarily and swap again after the last sqrt is computed.

                // First compute whether any value is zero
                let mut xi_list = [x0, x1, x2, x3, x4, x5, x6, x7];
                let mut zero_index = 0;
                let mut any_zero = 0;
                for (i, xi) in xi_list.iter().enumerate() {
                    // Keep track of whether any value is zero
                    let xi_is_zero = xi.iszero();
                    any_zero |= xi_is_zero;

                    // Compute the index which is zero
                    zero_index |= (i as u32 & xi_is_zero);
                }
                // If none of x0,...,x6 are zero, then set zero index to seven
                zero_index |= (7u32 & !any_zero);

                // Swap elements to ensure if there is a zero it's in seventh place
                xi_list.swap(zero_index as usize, 7);

                // Unpack the xi
                let [x0, x1, x2, x3, x4, x5, x6, x7] = xi_list;

                // Compute six square roots, y0 and y7 are set to ensure compatibility
                let y0;
                let mut y1 = (&x0 * &x1).sqrt().0;
                let mut y2 = (&x0 * &x2).sqrt().0;
                let mut y3 = (&x0 * &x3).sqrt().0;
                let mut y4 = (&x0 * &x4).sqrt().0;
                let mut y5 = (&x0 * &x5).sqrt().0;
                let mut y6 = (&x0 * &x6).sqrt().0;
                let y7;

                // Conditionally negate six square roots based on message bits
                let ctl1 = ((bits[0] as u32) & 1).wrapping_neg();
                y1.set_condneg(ctl1);

                let ctl2 = ((bits[1] as u32) & 1).wrapping_neg();
                y2.set_condneg(ctl2);

                let ctl3 = ((bits[2] as u32) & 1).wrapping_neg();
                y3.set_condneg(ctl3);

                let ctl4 = ((bits[3] as u32) & 1).wrapping_neg();
                y4.set_condneg(ctl4);

                let ctl5 = ((bits[4] as u32) & 1).wrapping_neg();
                y5.set_condneg(ctl5);

                let ctl6 = ((bits[5] as u32) & 1).wrapping_neg();
                y6.set_condneg(ctl6);

                // Compute the last square y7 given y0,...,y6
                (y0, y1, y2, y3, y4, y5, y6, y7) = self.last_sqrt(
                    &x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7, x0, y1, y2, y3, y4, y5, y6,
                );

                // Now we have computed the last square-root, we must swap
                // coordinates in the case that one value of xi was zero.
                let mut yi_list = [y0, y1, y2, y3, y4, y5, y6, y7];
                yi_list.swap(zero_index as usize, 7);
                let [y0, y1, y2, y3, y4, y5, y6, y7] = yi_list;

                let (b0, b1, b2, b3, b4, b5, b6, b7) =
                    self.to_hadamard(&y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7);

                ThetaPointDim3::new(&b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7)
            }

            fn last_sqrt(
                self,
                x0: &Fq,
                x1: &Fq,
                x2: &Fq,
                x3: &Fq,
                x4: &Fq,
                x5: &Fq,
                x6: &Fq,
                x7: &Fq,
                mut y0: Fq,
                mut y1: Fq,
                mut y2: Fq,
                mut y3: Fq,
                mut y4: Fq,
                mut y5: Fq,
                mut y6: Fq,
            ) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let (b0, b1, b2, b3, b4, b5, b6, b7) =
                    self.to_hadamard(x0, x1, x2, x3, x4, x5, x6, x7);

                let R1 = &b0 * &b1 * &b2 * &b3;
                let R3 = &b4 * &b5 * &b6 * &b7;

                let x04 = x0 * x4;
                let x15 = x1 * x5;
                let x26 = x2 * x6;
                let x37 = x3 * x7;
                let x0246 = &x04 * &x26;
                let x1357 = &x15 * &x37;

                let tmp = (&x04 - &x15 + &x26 - &x37).square();
                let t0 = &tmp.mul_small(16) - &(&x0246 + &x1357).mul_small(64);

                let r0 = (&R1 + &R3 - &t0);
                let r0_is_zero = r0.iszero();

                let y = y1 * y2 * y3 * y4 * y5 * y6;

                let mut t1 = r0.square() + (x0246 * x1357).mul_small(16384) - (R1 * R3).mul4();
                let mut t2 = r0.mul_small(256) * y;

                // When r0 == 0, t1 needs to be different
                let (a0, a1, a2, a3, a4, a5, a6, a7) = self.coords();
                let t1_prime = -(&a0 * &a1 * &a2 * &a3 * &a4 * &a5 * &a6 * &a7).mul_small(64);

                // If r0 is zero, modify the values of t1, t2
                t1.set_cond(&t1_prime, r0_is_zero);
                t2.set_cond(&y, r0_is_zero);

                // Scale y0...y6 by the denominator
                y0 *= t2;
                y1 *= t2;
                y2 *= t2;
                y3 *= t2;
                y4 *= t2;
                y5 *= t2;
                y6 *= t2;

                // Set y7 to the numerator with scaling factor
                let y7 = &t1 * x0 * x0.square();

                (y0, y1, y2, y3, y4, y5, y6, y7)
            }

            pub fn to_hash(self) -> (Fq, Fq, Fq, Fq, Fq, Fq, Fq) {
                let (X, Z, U, V, G, H, I, J) = (
                    self.X, self.Z, self.U, self.V, self.G, self.H, self.I, self.J,
                );
                let X_inv = X.invert();

                (
                    &Z * &X_inv,
                    &U * &X_inv,
                    &V * &X_inv,
                    &G * &X_inv,
                    &H * &X_inv,
                    &I * &X_inv,
                    &J * &X_inv,
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
