use super::Pattern;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;

pub struct Gradient {
    color_from: Color,
    color_to: Color,
    inv_tf: Matrix,
}

impl Gradient {
    pub fn new(color_from: Color, color_to: Color) -> Self {
        Self {
            color_from,
            color_to,
            inv_tf: Matrix::eye(4),
        }
    }

    pub fn with_transform(mut self, tf: Matrix) -> Self {
        self.inv_tf = tf.inverse();
        self
    }
}

impl Pattern for Gradient {
    fn color_at(&self, p: &Point) -> Color {
        self.color_from + (self.color_to - self.color_from) * (p.x - f64::floor(p.x))
    }

    fn inverse_transform(&self) -> &Matrix {
        &self.inv_tf
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::color;

    #[test]
    fn color_at() {
        let pattern = Gradient::new(color::WHITE, color::BLACK);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), color::WHITE);
        assert_eq!(
            pattern.color_at(&Point::new(0.25, 0., 0.)),
            Color::new(0.75, 0.75, 0.75)
        );
        assert_eq!(
            pattern.color_at(&Point::new(0.5, 0., 0.)),
            Color::new(0.5, 0.5, 0.5)
        );
        assert_eq!(
            pattern.color_at(&Point::new(0.75, 0., 0.)),
            Color::new(0.25, 0.25, 0.25)
        );
    }
}
