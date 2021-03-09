use std::cmp::PartialEq;

use crate::intersection::{Intersection, IntersectionList};
use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::vector::Vector;

#[cfg_attr(test, derive(Debug, PartialEq))]
pub struct Sphere {
    material: Material,

    transform: Matrix,
    inv_transform: Matrix,
    transpose_inv: Matrix,
    center: Point,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::eye(4),
            inv_transform: Matrix::eye(4),
            transpose_inv: Matrix::eye(4),
            center: Point::new(0., 0., 0.),
            material: Material::default(),
        }
    }

    pub fn glass() -> Self {
        let mut m = Material::default();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        Self {
            transform: Matrix::eye(4),
            inv_transform: Matrix::eye(4),
            transpose_inv: Matrix::eye(4),
            center: Point::new(0., 0., 0.),
            material: m,
        }
    }

    pub fn set_transform(mut self, t: Matrix) -> Self {
        self.inv_transform = t.inverse();
        self.transpose_inv = self.inv_transform.clone().transpose();
        self.transform = t;

        self
    }

    pub fn set_material(mut self, m: Material) -> Self {
        self.material = m;
        self
    }
}

impl Shape for Sphere {
    #[inline(always)]
    fn transform(&self) -> &Matrix {
        &self.transform
    }

    #[inline(always)]
    fn inverse_transform(&self) -> &Matrix {
        &self.inv_transform
    }

    #[inline(always)]
    fn transpose_inverse(&self) -> &Matrix {
        &self.transpose_inv
    }

    #[inline(always)]
    fn material(&self) -> &Material {
        &self.material
    }

    fn local_intersect(&self, ray: Ray) -> Vec<f64> {
        let sphere_to_ray = ray.origin - Point::new(0., 0., 0.);

        let a = Vector::dot(&ray.direction, &ray.direction);
        let b = 2. * Vector::dot(&ray.direction, &sphere_to_ray);
        let c = Vector::dot(&sphere_to_ray, &sphere_to_ray) - 1.;
        let discr = b * b - 4. * a * c;

        if discr < 0. {
            return Vec::new();
        }

        let mut t1 = (-b - f64::sqrt(discr)) / (2. * a);
        let mut t2 = (-b + f64::sqrt(discr)) / (2. * a);
        if t1 > t2 {
            std::mem::swap(&mut t1, &mut t2);
        }

        vec![t1, t2]
    }

    fn local_normal_at(&self, p: Point) -> Vector {
        p - self.center
    }
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use super::*;
    use crate::matrix::Axis;

    #[test]
    fn transform() {
        let t = Matrix::translation(2., 3., 4.);
        let s = Sphere::new().set_transform(t.clone());
        assert_eq!(*s.transform(), t.clone());
        assert_eq!(*s.inverse_transform(), t.inverse());
    }

    #[test]
    fn normal_at() {
        // the normal on a sphere at a point on the x axis
        let s = Sphere::new();
        let n = s.normal_at(Point::new(1., 0., 0.));
        assert_eq!(n, Vector::new(1., 0., 0.));

        // the normal on a sphere at a point on the y axis
        let s = Sphere::new();
        let n = s.normal_at(Point::new(0., 1., 0.));
        assert_eq!(n, Vector::new(0., 1., 0.));

        // the normal on a sphere at a point on the z axis
        let s = Sphere::new();
        let n = s.normal_at(Point::new(0., 0., 1.));
        assert_eq!(n, Vector::new(0., 0., 1.));

        // the normal on a sphere at a nonaxial point
        let s = Sphere::new();
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
        let s = Sphere::new();
        let n = s.normal_at(Point::new(
            f64::sqrt(3.) / 3.,
            f64::sqrt(3.) / 3.,
            f64::sqrt(3.) / 3.,
        ));
        assert_eq!(n, n.normalize());

        // computing the normal of a translated sphere
        let s = Sphere::new().set_transform(Matrix::translation(0., 1., 0.));
        let n = s.normal_at(Point::new(0., 1.70711, -0.70711));
        assert_eq!(n, Vector::new(0., 0.707106781, -0.707106781));

        // computing the normal of a transformed sphere
        let t = Matrix::scaling(1., 0.5, 1.) * Matrix::rotation(Axis::Z, PI / 5.);
        let s = Sphere::new().set_transform(t);
        let n = s.normal_at(Point::new(0., f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.));
        assert_eq!(n, Vector::new(0., 0.970142500, -0.242535625));
    }

    #[test]
    fn material() {
        let s = Sphere::new();
        assert_eq!(s.material, Material::default());

        let mut s = Sphere::new();
        let mut m = Material::default();

        m.ambient = 1.;
        s.material = Material::default();
        assert_eq!(s.material, Material::default());
    }

    #[test]
    fn glass() {
        let s = Sphere::glass();
        assert_eq!(s.transform, Matrix::eye(4));
        assert_eq!(s.material.transparency, 1.0);
        assert_eq!(s.material.refractive_index, 1.5);
    }
}
