use std::cmp::PartialEq;
use std::ops::{Add, Mul, Sub};

const TOL: f64 = 1e-9;
pub const WHITE: Color = Color {
    r: 1.,
    g: 1.,
    b: 1.,
};
pub const BLACK: Color = Color {
    r: 0.,
    g: 0.,
    b: 0.,
};
pub const DBROWN: Color = Color {
    r: 0.2,
    g: 0.2,
    b: 0.1,
};
pub const LBROWN: Color = Color {
    r: 0.3,
    g: 0.3,
    b: 0.15,
};

#[derive(Clone, Copy)]
#[cfg_attr(test, derive(Debug))]
pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Color {
    pub fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b }
    }
}

impl PartialEq<Color> for Color {
    fn eq(&self, rhs: &Color) -> bool {
        f64::abs(self.r - rhs.r) < TOL
            && f64::abs(self.g - rhs.g) < TOL
            && f64::abs(self.b - rhs.b) < TOL
    }
}

impl Add<Color> for Color {
    type Output = Self;

    fn add(self, rhs: Color) -> Self::Output {
        Self::Output {
            r: self.r + rhs.r,
            g: self.g + rhs.g,
            b: self.b + rhs.b,
        }
    }
}

impl Mul<Color> for Color {
    type Output = Self;

    fn mul(self, rhs: Color) -> Self::Output {
        Self::Output {
            r: self.r * rhs.r,
            g: self.g * rhs.g,
            b: self.b * rhs.b,
        }
    }
}

impl Mul<f64> for Color {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            r: self.r * rhs,
            g: self.g * rhs,
            b: self.b * rhs,
        }
    }
}

impl Mul<Color> for f64 {
    type Output = Color;

    fn mul(self, rhs: Color) -> Self::Output {
        rhs * self
    }
}

impl Sub<Color> for Color {
    type Output = Self;

    fn sub(self, rhs: Color) -> Self::Output {
        Self::Output {
            r: self.r - rhs.r,
            g: self.g - rhs.g,
            b: self.b - rhs.b,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new() {
        let c = Color::new(-0.5, 0.4, 1.7);
        assert_eq!(c.r, -0.5);
        assert_eq!(c.g, 0.4);
        assert_eq!(c.b, 1.7);
    }

    #[test]
    fn add() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 + c2, Color::new(1.6, 0.7, 1.0));
    }

    #[test]
    fn mul() {
        let c = Color::new(0.2, 0.3, 0.4);
        assert_eq!(c * 2., Color::new(0.4, 0.6, 0.8));

        let c1 = Color::new(1., 0.2, 0.4);
        let c2 = Color::new(0.9, 1., 0.1);
        assert_eq!(c1 * c2, Color::new(0.9, 0.2, 0.04));
    }

    #[test]
    fn sub() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 - c2, Color::new(0.2, 0.5, 0.5));
    }
}
