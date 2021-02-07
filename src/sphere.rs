use std::cmp::PartialEq;

use crate::matrix::Matrix;

#[derive(Clone, Debug)]
pub struct Sphere {
    pub transform: Matrix,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::eye(4),
        }
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
        s.transform = t.clone();
        assert_eq!(s.transform, t);
    }
}
