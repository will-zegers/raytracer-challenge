use std::cmp::PartialEq;

use crate::matrix::Matrix;

#[derive(Clone, Debug)]
pub struct Sphere {
    transform: Matrix,
    inv_transform: Matrix,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::eye(4),
            inv_transform: Matrix::eye(4),
        }
    }

    pub fn set_transform(&mut self, t: Matrix) {
        self.inv_transform = t.inverse();
        self.transform = t;
    }

    pub fn transform(&self) -> &Matrix {
        &self.transform
    }

    pub fn inverse_transform(&self) -> &Matrix {
        &self.inv_transform
    }
}

impl PartialEq<Sphere> for Sphere {
    fn eq(&self, _rhs: &Self) -> bool {
        true
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn transform() {
        let mut s = Sphere::new();
        let t = Matrix::translation(2., 3., 4.);
        s.set_transform(t.clone());
        assert_eq!(*s.transform(), t.clone());
        assert_eq!(*s.inverse_transform(), t.inverse());
    }
}
