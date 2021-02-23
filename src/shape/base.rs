use crate::matrix::Matrix;

pub trait Shape {
    fn transform(&self) -> &Matrix;
    fn inverse_transform(&self) -> &Matrix;
}
