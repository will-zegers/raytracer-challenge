use super::Shape;

use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

pub struct Plane {
    material: Material,

    transform: Matrix,
    inv_transform: Matrix,
    transpose_inv: Matrix,
}

impl Plane {
    pub fn new(transform: Matrix, material: Material) -> Self {
        let inv_transform = transform.inverse();
        Self {
            transform,
            inv_transform: inv_transform.clone(),
            transpose_inv: inv_transform.transpose(),
            material,
        }
    }

    #[cfg(test)]
    pub fn default() -> Self {
        Self {
            transform: Matrix::eye(4),
            inv_transform: Matrix::eye(4),
            transpose_inv: Matrix::eye(4),
            material: Material::default(),
        }
    }

    pub fn set_material(mut self, m: Material) -> Self {
        self.material = m;
        self
    }

    pub fn set_transform(mut self, t: Matrix) -> Self {
        let inv_transform = t.inverse();

        self.transform = t;
        self.inv_transform = inv_transform.clone();
        self.transpose_inv = inv_transform.transpose();

        self
    }
}

impl Shape for Plane {
    fn transform(&self) -> &Matrix {
        &self.transform
    }

    fn inverse_transform(&self) -> &Matrix {
        &self.inv_transform
    }

    fn transpose_inverse(&self) -> &Matrix {
        &self.transpose_inv
    }

    fn material(&self) -> &Material {
        &self.material
    }

    fn local_normal_at(&self, _: Point) -> Vector {
        Vector::new(0., 1., 0.)
    }

    fn local_intersect(&self, r: Ray) -> Vec<f64> {
        const TOL: f64 = 1e-9;
        if f64::abs(r.direction.y) < TOL {
            return Vec::new();
        }

        let t = -r.origin.y / r.direction.y;
        vec![t]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn normal() {
        let p = Plane::default();
        let n1 = p.local_normal_at(Point::new(0., 0., 0.));
        let n2 = p.local_normal_at(Point::new(10., 0., -10.));
        let n3 = p.local_normal_at(Point::new(-5., 0., 150.));
        assert_eq!(n1, Vector::new(0., 1., 0.));
        assert_eq!(n2, Vector::new(0., 1., 0.));
        assert_eq!(n3, Vector::new(0., 1., 0.));
    }

    #[test]
    fn local_intersect() {
        // Intersect with a ray parallel to the line
        let p = Plane::default();
        let r = Ray::new(Point::new(0., 10., 0.), Vector::new(0., 0., 1.));
        let xs = p.local_intersect(r);
        assert!(xs.is_empty());

        // Intersect with a ray coplanar to the line
        let p = Plane::default();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let xs = p.local_intersect(r);
        assert!(xs.is_empty());

        // A ray intersecting a plane from above
        let p = Plane::default();
        let r = Ray::new(Point::new(0., 1., 0.), Vector::new(0., -1., 0.));
        let xs = p.local_intersect(r);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0], 1.);

        // A ray intersecting a plane from above
        let p = Plane::default();
        let r = Ray::new(Point::new(0., -1., 0.), Vector::new(0., 1., 0.));
        let xs = p.local_intersect(r);
        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0], 1.);
    }
}
