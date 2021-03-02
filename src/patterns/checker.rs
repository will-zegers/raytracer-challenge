use super::Pattern;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;

pub struct Checker {
    color1: Color,
    color2: Color,
    inv_tf: Matrix,
}

impl Checker {
    pub fn new(color1: Color, color2: Color) -> Self {
        Self {
            color1,
            color2,
            inv_tf: Matrix::eye(4),
        }
    }

    pub fn set_transform(mut self, tf: Matrix) -> Self {
        self.inv_tf = tf.inverse();

        self
    }
}

impl Pattern for Checker {
    fn color_at(&self, p: &Point) -> Color {
        let cumul_floor = f64::floor(p.x) + f64::floor(p.y) + f64::floor(p.z);
        if cumul_floor % 2. == 0. {
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
    use std::rc::Rc;

    use super::*;

    use crate::color::{BLACK, WHITE};

    #[test]
    fn color_at() {
        let pattern = Checker::new(WHITE, BLACK);

        // checkers should repeat in x
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0.99, 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(1.01, 0., 0.)), BLACK);

        // checkers should repeat in y
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0.99, 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 1.01, 0.)), BLACK);

        // checkers should repeat in z
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.99)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 1.01)), BLACK);
    }
}
