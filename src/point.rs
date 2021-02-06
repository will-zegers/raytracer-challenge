use std::cmp::PartialEq;
use std::ops::{Add, Div, Index, Mul, Neg, Sub};

use crate::vector::Vector;

const TOL: f64 = 1e-9;

#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    w: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 1. }
    }
}

impl PartialEq<Point> for Point {
    fn eq(&self, rhs: &Point) -> bool {
        self.x - rhs.x < TOL && self.y - rhs.y < TOL && self.z - rhs.z < TOL
    }
}

impl Add<Vector> for Point {
    type Output = Self;

    fn add(self, rhs: Vector) -> Self::Output {
        Self::Output {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: 1.,
        }
    }
}

impl Div<f64> for Point {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        (1. / rhs) * self
    }
}

impl Index<usize> for Point {
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

impl Mul<f64> for Point {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: rhs * self.x,
            y: rhs * self.y,
            z: rhs * self.z,
            w: 1.,
        }
    }
}

impl Mul<Point> for f64 {
    type Output = Point;

    fn mul(self, rhs: Point) -> Self::Output {
        rhs * self
    }
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Output {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: 1.,
        }
    }
}

impl Sub<Point> for Point {
    type Output = Vector;

    fn sub(self, rhs: Point) -> Self::Output {
        Self::Output::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Sub<Vector> for Point {
    type Output = Point;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: 1.,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn add() {
        let a1 = Point::new(3., -2., 5.);
        let a2 = Vector::new(-2., 3., 1.);
        assert_eq!(a1 + a2, Point::new(1., 1., 6.));
    }

    #[test]
    fn sub() {
        let p1 = Point::new(3., 2., 1.);
        let p2 = Point::new(5., 6., 7.);
        assert_eq!(p1 - p2, Vector::new(-2., -4., -6.));

        let p = Point::new(3., 2., 1.);
        let v = Vector::new(5., 6., 7.);
        assert_eq!(p - v, Point::new(-2., -4., -6.));
    }

    #[test]
    fn neg() {
        let a = Point::new(1., -2., 3.);
        assert_eq!(-a, Point::new(-1., 2., -3.));
    }

    #[test]
    fn mul() {
        let a = Point::new(1., -2., 3.);
        assert_eq!(a * 3.5, Point::new(3.5, -7., 10.5));

        let a = Point::new(1., -2., 3.);
        assert_eq!(a * 0.5, Point::new(0.5, 1., 1.5));
    }
}
