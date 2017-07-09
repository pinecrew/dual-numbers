use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::ops::{Add, Sub, Mul, Div, Neg};
use std::f64;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual<T> {
    real: f64,
    dual: T,
}

impl<T> Dual<T> {
    pub fn new(real: f64, dual: T) -> Dual<T> {
        Dual { real, dual }
    }

    pub fn sin(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        Dual::new(self.real.sin(), self.dual * self.real.cos())
    }

    pub fn cos(self) -> Dual<T>
    where
        T: Mul<f64, Output = T> + Neg<Output = T>,
    {
        Dual::new(self.real.cos(), -self.dual * self.real.sin())
    }

    pub fn exp(self) -> Dual<T>
    where
        T: Mul<f64, Output = T>,
    {
        let real_exp = self.real.exp();
        Dual::new(real_exp, self.dual * real_exp)
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
