#![cfg(test)]
use super::Pattern;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;

pub struct TestPattern {
    inv_tf: Matrix,
}

impl TestPattern {
    pub fn new() -> Self {
        Self {
            inv_tf: Matrix::eye(4),
        }
    }

    pub fn set_transform(mut self, tf: Matrix) -> Self {
        self.inv_tf = tf.inverse();
        self
    }
}

impl Pattern for TestPattern {
    fn color_at(&self, p: &Point) -> Color {
        Color::new(p.x, p.y, p.z)
    }

    fn inverse_transform(&self) -> &Matrix {
        &self.inv_tf
    }
}
