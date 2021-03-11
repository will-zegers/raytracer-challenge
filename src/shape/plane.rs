use super::Shape;

use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

pub struct Plane {
    material: Material,

    inv_tf: Matrix,
    xpose_inv_tf: Matrix,
}

impl Plane {
    pub fn new(transform: Matrix, material: Material) -> Self {
        let inv_tf = transform.inverse();
        Self {
            inv_tf: inv_tf.clone(),
            xpose_inv_tf: inv_tf.transpose(),
            material,
        }
    }

    #[cfg(test)]
    pub fn default() -> Self {
        Self {
            inv_tf: Matrix::eye(4),
            xpose_inv_tf: Matrix::eye(4),
            material: Material::default(),
        }
    }

    pub fn set_material(mut self, m: Material) -> Self {
        self.material = m;
        self
    }

    pub fn set_transform(mut self, t: Matrix) -> Self {
        let inv_tf = t.inverse();

        self.inv_tf = inv_tf.clone();
        self.xpose_inv_tf = inv_tf.transpose();

        self
    }
}

impl Shape for Plane {
    fn inverse_transform(&self) -> &Matrix {
        &self.inv_tf
    }

    fn transpose_inverse(&self) -> &Matrix {
        &self.xpose_inv_tf
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
