use std::ops::{Add, Div, Index, Mul, Neg, Sub};

use crate::point::Point;

#[derive(Clone, Copy)]
#[cfg_attr(test, derive(Debug))]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    w: f64,
}

// statics
impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 0. }
    }

    pub fn cross(v1: &Vector, v2: &Vector) -> Self {
        Self {
            x: v1.y * v2.z - v1.z * v2.y,
            y: v1.z * v2.x - v1.x * v2.z,
            z: v1.x * v2.y - v1.y * v2.x,
            w: 0.,
        }
    }

    pub fn dot(v1: &Vector, v2: &Vector) -> f64 {
        v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
    }
}

// methods
impl Vector {
    pub fn length(&self) -> f64 {
        f64::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }

    pub fn normalize(self) -> Self {
        (1. / self.length()) * self
    }

    pub fn reflect(self, normal: &Self) -> Self {
        self - 2. * Vector::dot(&self, normal) * *normal
    }
}

impl Add<Point> for Vector {
    type Output = Point;

    fn add(self, rhs: Point) -> Self::Output {
        rhs + self
    }
}

impl Div<f64> for Vector {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        (1. / rhs) * self
    }
}

impl Index<usize> for Vector {
    type Output = f64;

    fn index(&self, i: usize) -> &Self::Output {
        match i {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            3 => &self.w,
            _ => panic!("invalid vector index: {:?}", i),
        }
    }
}

impl Mul<f64> for Vector {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: rhs * self.x,
            y: rhs * self.y,
            z: rhs * self.z,
            w: self.w,
        }
    }
}

impl Mul<Vector> for f64 {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        rhs * self
    }
}

impl Neg for Vector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Output {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }
}

impl Sub<Vector> for Vector {
    type Output = Vector;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const TOL: f64 = 1e-9;

    impl PartialEq<Vector> for Vector {
        fn eq(&self, rhs: &Vector) -> bool {
            f64::abs(self.x - rhs.x) < TOL
                && f64::abs(self.y - rhs.y) < TOL
                && f64::abs(self.z - rhs.z) < TOL
        }
    }

    #[test]
    fn add() {
        let a1 = Vector::new(3., -2., 5.);
        let a2 = Point::new(-2., 3., 1.);
        assert_eq!(a1 + a2, Point::new(1., 1., 6.));
    }

    #[test]
    fn sub() {
        let v1 = Vector::new(3., 2., 1.);
        let v2 = Vector::new(5., 6., 7.);
        assert_eq!(v1 - v2, Vector::new(-2., -4., -6.));

        let v1 = Vector::new(0., 0., 0.);
        let v2 = Vector::new(1., -2., 3.);
        assert_eq!(v1 - v2, Vector::new(-1., 2., -3.));
    }

    #[test]
    fn neg() {
        let a = Vector::new(1., -2., 3.);
        assert_eq!(-a, Vector::new(-1., 2., -3.));
    }

    #[test]
    fn mul() {
        let a = Vector::new(1., -2., 3.);
        assert_eq!(a * 3.5, Vector::new(3.5, -7., 10.5));

        let a = Vector::new(1., -2., 3.);
        assert_eq!(a * 0.5, Vector::new(0.5, -1., 1.5));
    }

    #[test]
    fn div() {
        let a = Vector::new(1., -2., 3.);
        assert_eq!(a / 2., Vector::new(0.5, -1., 1.5));
    }

    #[test]
    fn length() {
        let v = Vector::new(1., 0., 0.);
        assert_eq!(v.length(), 1.);

        let v = Vector::new(0., 1., 0.);
        assert_eq!(v.length(), 1.);

        let v = Vector::new(0., 0., 1.);
        assert_eq!(v.length(), 1.);

        let v = Vector::new(1., 2., 3.);
        assert_eq!(v.length(), f64::sqrt(14.));

        let v = Vector::new(-1., -2., -3.);
        assert_eq!(v.length(), f64::sqrt(14.));
    }

    #[test]
    fn normalize() {
        let v = Vector::new(4., 0., 0.);
        assert_eq!(v.normalize(), Vector::new(1., 0., 0.));

        let v = Vector::new(1., 2., 3.);
        assert_eq!(
            v.normalize(),
            Vector::new(0.267261241, 0.534522483, 0.801783725)
        );

        let v = Vector::new(1., 2., 3.);
        assert_eq!(v.normalize().length(), 1.,);
    }

    #[test]
    fn dot() {
        let a = Vector::new(1., 2., 3.);
        let b = Vector::new(2., 3., 4.);
        assert_eq!(Vector::dot(&a, &b), 20.);
    }

    #[test]
    fn cross() {
        let a = Vector::new(1., 2., 3.);
        let b = Vector::new(2., 3., 4.);
        assert_eq!(Vector::cross(&a, &b), Vector::new(-1., 2., -1.));
        assert_eq!(Vector::cross(&b, &a), Vector::new(1., -2., 1.));
    }

    #[test]
    fn reflect() {
        // reflecting a vector approaching at 45⁰
        let v = Vector::new(1., -1., 0.);
        let n = Vector::new(0., 1., 0.);
        let r = v.reflect(&n);
        assert_eq!(r, Vector::new(1., 1., 0.));

        // reflecting a vector off a slanted surface
        let v = Vector::new(0., -1., 0.);
        let n = Vector::new(f64::sqrt(2.) / 2., f64::sqrt(2.) / 2., 0.);
        let r = v.reflect(&n);
        assert_eq!(r, Vector::new(1., 0., 0.));
    }
}
