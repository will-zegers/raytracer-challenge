use super::Pattern;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;

pub struct Stripe {
    color1: Color,
    color2: Color,
    inv_tf: Matrix,
}

impl Stripe {
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

impl Pattern for Stripe {
    fn color_at(&self, p: &Point) -> Color {
        if f64::floor(p.x) % 2. == 0. {
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
    use crate::material::Material;
    use crate::shape::Sphere;

    #[test]
    fn new() {
        let pattern = Stripe::new(WHITE, BLACK);
        assert_eq!(pattern.color1, WHITE);
        assert_eq!(pattern.color2, BLACK);
    }

    #[test]
    fn color_at() {
        let pattern = Stripe::new(WHITE, BLACK);

        // a stripe pattern is constant in y
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 1., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 2., 0.)), WHITE);

        // a stripe pattern is constant in z
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 1.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0., 0., 2.)), WHITE);

        // a stripe pattern alternates in x
        assert_eq!(pattern.color_at(&Point::new(0., 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(0.9, 0., 0.)), WHITE);
        assert_eq!(pattern.color_at(&Point::new(1., 0., 0.)), BLACK);
        assert_eq!(pattern.color_at(&Point::new(-0.1, 0., 0.)), BLACK);
        assert_eq!(pattern.color_at(&Point::new(-1., 0., 0.)), BLACK);
    }

    #[test]
    fn stripe_at_object() {
        // stripes with object transformation
        let object = Sphere::new().with_transform(Matrix::scaling(2., 2., 2.));
        let object = Rc::new(object);
        let pattern: Box<dyn Pattern> = Box::new(Stripe::new(WHITE, BLACK));
        assert_eq!(
            pattern.color_at_object(object, &Point::new(1.5, 0., 0.)),
            WHITE
        );

        // stripes with pattern transformation
        let object = Sphere::new();
        let object = Rc::new(object);
        let tf = Matrix::scaling(2., 2., 2.);
        let pattern: Box<dyn Pattern> = Box::new(Stripe::new(WHITE, BLACK).with_transform(tf));
        assert_eq!(
            pattern.color_at_object(object, &Point::new(1.5, 0., 0.)),
            WHITE
        );

        let object = Sphere::new().with_transform(Matrix::scaling(2., 2., 2.));
        let object = Rc::new(object);
        let pattern: Box<dyn Pattern> =
            Box::new(Stripe::new(WHITE, BLACK).with_transform(Matrix::translation(0.5, 0., 0.)));
        assert_eq!(
            pattern.color_at_object(object, &Point::new(3.5, 0., 0.)),
            BLACK
        );
    }
}
