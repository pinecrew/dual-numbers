use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::ops::{Add, Sub, Mul, Div, Neg};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual<T> {
    real: f64,
    dual: T,
}

impl<T> Dual<T> {
    pub fn new<U: Into<f64>>(real: U, dual: T) -> Dual<T> {
        Dual {
            real: real.into(),
            dual: dual,
        }
    }

    pub fn sin(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(self.real.sin(), self.dual * self.real.cos())
    }

    pub fn cos(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(self.real.cos(), self.dual * (-self.real.sin()))
    }

    pub fn tan(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.tan(), self.dual / self.real.cos().powi(2))
    }

    pub fn asin(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(
            self.real.asin(),
            self.dual / (1f64 - self.real.powi(2)).sqrt(),
        )
    }

    pub fn acos(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(
            self.real.acos(),
            self.dual / (-(1f64 - self.real.powi(2)).sqrt()),
        )
    }

    pub fn atan(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.atan(), self.dual / (1f64 + self.real.powi(2)))
    }

    pub fn sinh(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(self.real.sinh(), self.dual * self.real.cosh())
    }

    pub fn cosh(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(self.real.cosh(), self.dual * self.real.sinh())
    }

    pub fn tanh(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.tanh(), self.dual / self.real.cosh().powi(2))
    }

    pub fn asinh(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(
            self.real.asinh(),
            self.dual / (self.real.powi(2) + 1f64).sqrt(),
        )
    }

    pub fn acosh(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(
            self.real.acosh(),
            self.dual / (self.real.powi(2) - 1f64).sqrt(),
        )
    }

    pub fn atanh(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.atanh(), self.dual / (1f64 - self.real.powi(2)))
    }

    pub fn exp(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        let real_exp = self.real.exp();
        Dual::new(real_exp, self.dual * real_exp)
    }

    pub fn exp2(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        let real_exp2 = self.real.exp2();
        Dual::new(real_exp2, self.dual * real_exp2 * 2f64.ln())
    }

    pub fn ln(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.ln(), self.dual / self.real)
    }

    pub fn log(self, base: f64) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.log(base), self.dual / self.real / base.ln())
    }

    pub fn log2(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.log2(), self.dual / self.real / 2f64.ln())
    }

    pub fn log10(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.log10(), self.dual / self.real / 10f64.ln())
    }

    pub fn abs(self) -> f64 {
        self.real
    }

    pub fn powi(self, n: i32) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(
            self.real.powi(n),
            self.dual * (self.real.powi(n - 1) * n as f64),
        )
    }

    pub fn powf(self, n: f64) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(
            self.real.powf(n),
            self.dual * (self.real.powf(n - 1f64) * n),
        )
    }

    pub fn sqrt(self) -> Dual<T>
    where
        T: Div<f64, Output = T>,
    {
        Dual::new(self.real.sqrt(), self.dual / self.real.sqrt() / 2f64)
    }
}

impl<T> Neg for Dual<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Dual::new(-self.real, -self.dual)
    }
}

impl<T> Add for Dual<T>
where
    T: AddAssign,
{
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T> AddAssign for Dual<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, other: Self) {
        self.real += other.real;
        self.dual += other.dual;
    }
}

impl<T> Sub for Dual<T>
where
    T: SubAssign,
{
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T> SubAssign for Dual<T>
where
    T: SubAssign,
{
    fn sub_assign(&mut self, other: Self) {
        self.real -= other.real;
        self.dual -= other.dual;
    }
}

impl<T> Mul for Dual<T>
where
    T: Mul<f64, Output = T> + Add<T, Output = T> + Copy,
{
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<T, U> Mul<U> for Dual<T>
where
    T: MulAssign<f64>,
    U: Into<f64>,
{
    type Output = Self;

    fn mul(mut self, rhs: U) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<T> MulAssign for Dual<T>
where
    T: Mul<f64, Output = T> + Add<T, Output = T> + Copy,
{
    fn mul_assign(&mut self, other: Self) {
        self.dual = other.dual * self.real + self.dual * other.real;
        self.real *= other.real;
    }
}

impl<T, U> MulAssign<U> for Dual<T>
where
    T: MulAssign<f64>,
    U: Into<f64>,
{
    fn mul_assign(&mut self, other: U) {
        let tmp: f64 = other.into();
        self.real *= tmp;
        self.dual *= tmp;
    }
}

impl<T> Div for Dual<T>
where
    T: Mul<f64, Output = T> + Div<f64, Output = T> + Sub<T, Output = T> + Copy,
{
    type Output = Self;

    fn div(mut self, rhs: Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<T, U> Div<U> for Dual<T>
where
    T: DivAssign<f64>,
    U: Into<f64>,
{
    type Output = Self;

    fn div(mut self, rhs: U) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<T> DivAssign for Dual<T>
where
    T: Mul<f64, Output=T> + Div<f64, Output=T> + Sub<T, Output=T> + Copy
{
    fn div_assign(&mut self, other: Self) {
        self.dual = other.dual / self.real - self.dual * other.real / self.real.powi(2);
        self.real /= other.real;
    }
}

impl<T, U> DivAssign<U> for Dual<T>
where
    T: DivAssign<f64>,
    U: Into<f64>,
{
    fn div_assign(&mut self, other: U) {
        let tmp: f64 = other.into();
        self.real /= tmp;
        self.dual /= tmp;
    }
}

#[cfg(test)]
mod dual_test {
    use super::*;
    use std::f64::consts::*;

    // change this
    const EPS: f64 = 1E-15;

    #[test]
    fn add() {
        let a = Dual::new(1_f64, 2_f64);
        let b = Dual::new(3_f64, 4_f64);
        let c = Dual::new(4_f64, 6_f64);
        assert_eq!(a + b, c);
    }

    #[test]
    fn sub() {
        let a = Dual::new(1_f64, 2_f64);
        let b = Dual::new(3_f64, 4_f64);
        let c = Dual::new(-2_f64, -2_f64);
        assert_eq!(a - b, c);
    }

    #[test]
    fn mul() {
        let a = Dual::new(1_f64, 4_f64);
        let b = Dual::new(5_f64, 2_f64);
        let c = Dual::new(5_f64, 22_f64);
        assert_eq!(a * b, c);
    }

    #[test]
    fn div() {
        let a = Dual::new(4_f64, 3_f64);
        let b = Dual::new(1_f64, 2_f64);
        let c = Dual::new(4_f64, 0.3125_f64);
        assert_eq!(a / b, c);
    }

    #[test]
    fn neg() {
        let a = Dual::new(-1_f64, 1_f64);
        let b = Dual::new(1_f64, -1_f64);
        assert_eq!(-a, b);
    }

    #[test]
    fn sin() {
        let a = Dual::new(FRAC_PI_2, 2_f64);
        let b = Dual::new(1_f64, 0_f64);
        let c = a.sin() - b;
        assert_eq!(c.real == 0_f64 && c.dual.abs() < EPS, true);
    }

    #[test]
    fn cos() {
        let a = Dual::new(2_f64 * PI, 2_f64);
        let b = Dual::new(1_f64, 0_f64);
        let c = a.cos() - b;
        assert_eq!(c.real == 0_f64 && c.dual.abs() < EPS, true);
    }

    #[test]
    fn tan() {
        let a = Dual::new(FRAC_PI_4, 0.5_f64);
        let b = Dual::new(1_f64, 1_f64);
        let c = a.tan() - b;
        assert_eq!(c.real.abs() < EPS && c.dual.abs() < EPS, true);
    }

    // #[test]
    // fn asin() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn acos() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn atan() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn sinh() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn cosh() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn tanh() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn asinh() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn acosh() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn atanh() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    #[test]
    fn exp() {
        let a = Dual::new(0_f64, 2_f64);
        let b = Dual::new(1_f64, 2_f64);
        assert_eq!(a.exp(), b);
    }

    #[test]
    fn exp2() {
        let a = Dual::new(0_f64, 2_f64);
        let b = Dual::new(1_f64, 2_f64 * 2_f64.ln());
        assert_eq!(a.exp2(), b);
    }

    // #[test]
    // fn ln() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn log() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn log2() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn log10() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn abs() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn powi() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn powf() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }

    // #[test]
    // fn sqrt() {
    //     let a = Dual::new(0_f64, 0_f64);
    //     let b = Dual::new(0_f64, 0_f64);
    //     assert_eq!(true, true);
    // }
}
