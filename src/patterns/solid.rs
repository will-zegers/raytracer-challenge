use super::Pattern;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;

pub struct Solid {
    color: Color,
    inv_tf: Matrix,
}

impl Solid {
    pub fn new(color: Color) -> Self {
        Self {
            color,
            inv_tf: Matrix::eye(4),
        }
    }

    pub fn with_transform(mut self, tf: Matrix) -> Self {
        self.inv_tf = tf.inverse();
        self
    }
}

impl Pattern for Solid {
    #[inline(always)]
    fn color_at(&self, _: &Point) -> Color {
        self.color
    }

    fn inverse_transform(&self) -> &Matrix {
        &self.inv_tf
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::color::WHITE;

    #[test]
    fn new() {
        let pattern = Solid::new(WHITE);
        assert_eq!(pattern.color, WHITE);
    }

    #[test]
    fn color_at() {
        let pattern = Solid::new(WHITE);

        // a solid pattern is constant in y
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 1., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 2., 0.)), WHITE);

        // a solid pattern is constant in z
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 1.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 2.)), WHITE);

        // a solid pattern is constant in x
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0.9, 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(1., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(-0.1, 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(-1., 0., 0.)), WHITE);
    }
}
