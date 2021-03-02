use std::fmt::Debug;
use std::rc::Rc;

use crate::color::Color;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;

pub trait Pattern {
    fn color_at(&self, _: &Point) -> Color;
    fn inverse_transform(&self) -> &Matrix;

    fn color_at_object(&self, object: Rc<dyn Shape>, p: &Point) -> Color {
        let object_p = object.inverse_transform() * *p;
        let pattern_p = self.inverse_transform() * object_p;
        self.color_at(&pattern_p)
    }
}

#[cfg(test)]
impl Debug for dyn Pattern {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        Ok(())
    }
}
