use super::Pattern;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;

pub struct Ring {
    color1: Color,
    color2: Color,
    inv_tf: Matrix,
}

impl Ring {
    pub fn new(color1: Color, color2: Color) -> Self {
        Self {
            color1,
            color2,
            inv_tf: Matrix::eye(4),
        }
    }

    pub fn with_transform(mut self, tf: Matrix) -> Self {
        self.inv_tf = tf.inverse();
        self
    }
}

impl Pattern for Ring {
    fn color_at(&self, p: &Point) -> Color {
        let r = f64::sqrt(p.x * p.x + p.z * p.z);
        if f64::floor(r) % 2. == 0. {
            self.color1
        } else {
            self.color2
        }
    }

    fn inverse_transform(&self) -> &Matrix {
        &self.inv_tf
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::color::{BLACK, WHITE};

    #[test]
    fn color_at() {
        let pattern = Ring::new(WHITE, BLACK);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(1., 0., 0.)), BLACK);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 1.)), BLACK);
        assert_eq!(
            pattern.color_at(&Point::new(f64::sqrt(2.) / 2., 0., f64::sqrt(2.) / 2.)),
            BLACK
        );
    }
}
