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

impl<T> MulAssign for Dual<T>
where
    T: Mul<f64, Output = T> + Add<T, Output = T> + Copy,
{
    fn mul_assign(&mut self, other: Self) {
        self.dual = other.dual * self.real + self.dual * other.real;
        self.real *= other.real;
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

impl<T> DivAssign for Dual<T>
where
    T: Mul<f64, Output=T> + Div<f64, Output=T> + Sub<T, Output=T> + Copy
{
    fn div_assign(&mut self, other: Self) {
        self.dual = other.dual / self.real - self.dual * other.real / self.real.powi(2);
        self.real /= other.real;
    }
}

#[cfg(test)]
mod dual_test {
    use super::*;


    #[test]
    fn mul() {
        let a = Dual::new(1f64, 4f64);
        let b = Dual::new(5f64, 2f64);
        assert_eq!(a * b, Dual::new(5f64, 22f64));
    }
}
