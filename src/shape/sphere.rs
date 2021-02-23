use std::cmp::PartialEq;

use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;
use crate::vector::Vector;

#[derive(Clone, Debug, PartialEq)]
pub struct Sphere {
    pub material: Material,

    transform: Matrix,
    inv_transform: Matrix,
    transpose_inv: Matrix,
    center: Point,
}

impl Sphere {
    pub fn new(transform: Matrix, material: Material) -> Self {
        let inv_transform = transform.inverse();
        let transpose_inv = inv_transform.clone().transpose();
        Self {
            transform,
            inv_transform,
            transpose_inv,
            center: Point::new(0., 0., 0.),
            material,
        }
    }

    pub fn default() -> Self {
        Self {
            transform: Matrix::eye(4),
            inv_transform: Matrix::eye(4),
            transpose_inv: Matrix::eye(4),
            center: Point::new(0., 0., 0.),
            material: Material::default(),
        }
    }

    pub fn set_transform(&mut self, t: Matrix) {
        self.inv_transform = t.inverse();
        self.transpose_inv = self.inv_transform.clone().transpose();
        self.transform = t;
    }

    pub fn normal_at(&self, p: Point) -> Vector {
        let obj_point = &self.inv_transform * p;
        let obj_normal = obj_point - self.center;
        let world_normal = &self.transpose_inv * obj_normal;
        world_normal.normalize()
    }
}

impl Shape for Sphere {
    #[inline(always)]
    fn transform(&self) -> &Matrix {
        &self.transform
    }

    #[allow(dead_code)]
    fn inverse_transform(&self) -> &Matrix {
        &self.inv_transform
    }
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use super::*;
    use crate::matrix::Axis;

    #[test]
    fn transform() {
        let mut s = Sphere::default();
        let t = Matrix::translation(2., 3., 4.);
        s.set_transform(t.clone());
        assert_eq!(*s.transform(), t.clone());
        assert_eq!(*s.inverse_transform(), t.inverse());
    }

    #[test]
    fn normal_at() {
        // the normal on a sphere at a point on the x axis
        let s = Sphere::default();
        let n = s.normal_at(Point::new(1., 0., 0.));
        assert_eq!(n, Vector::new(1., 0., 0.));

        // the normal on a sphere at a point on the y axis
        let s = Sphere::default();
        let n = s.normal_at(Point::new(0., 1., 0.));
        assert_eq!(n, Vector::new(0., 1., 0.));

        // the normal on a sphere at a point on the z axis
        let s = Sphere::default();
        let n = s.normal_at(Point::new(0., 0., 1.));
        assert_eq!(n, Vector::new(0., 0., 1.));

        // the normal on a sphere at a nonaxial point
        let s = Sphere::default();
        let n = s.normal_at(Point::new(
            f64::sqrt(3.) / 3.,
            f64::sqrt(3.) / 3.,
            f64::sqrt(3.) / 3.,
        ));
        assert_eq!(
            n,
            Vector::new(f64::sqrt(3.) / 3., f64::sqrt(3.) / 3., f64::sqrt(3.) / 3.)
        );

        // normals are all normalized
        let s = Sphere::default();
        let n = s.normal_at(Point::new(
            f64::sqrt(3.) / 3.,
            f64::sqrt(3.) / 3.,
            f64::sqrt(3.) / 3.,
        ));
        assert_eq!(n, n.normalize());

        // computing the normal of a translated sphere
        let mut s = Sphere::default();
        s.set_transform(Matrix::translation(0., 1., 0.));
        let n = s.normal_at(Point::new(0., 1.70711, -0.70711));
        assert_eq!(n, Vector::new(0., 0.707106781, -0.707106781));

        // computing the normal of a transformed sphere
        let mut s = Sphere::default();
        let t = Matrix::scaling(1., 0.5, 1.) * Matrix::rotation(Axis::Z, PI / 5.);
        s.set_transform(t);
        let n = s.normal_at(Point::new(0., f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.));
        assert_eq!(n, Vector::new(0., 0.970142500, -0.242535625));
    }

    #[test]
    fn material() {
        let s = Sphere::default();
        assert_eq!(s.material, Material::default());

        let mut s = Sphere::default();
        let mut m = Material::default();

        m.ambient = 1.;
        s.material = m.clone();
        assert_eq!(s.material, m);
    }
}
