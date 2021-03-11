use std::fmt::Debug;

use crate::intersection::IntersectionList;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

pub trait Shape {
    fn inverse_transform(&self) -> &Matrix;
    fn transpose_inverse(&self) -> &Matrix;
    fn material(&self) -> &Material;

    fn local_intersect(&self, _: Ray) -> Vec<f64>;
    fn local_normal_at(&self, _: Point) -> Vector;

    fn normal_at(&self, p: Point) -> Vector {
        let local_point = self.inverse_transform() * p;
        let local_normal = self.local_normal_at(local_point);
        let world_normal = self.transpose_inverse() * local_normal;

        world_normal.normalize()
    }
}

#[cfg(test)]
impl Debug for dyn Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self.material())
    }
}
