use std::f64::EPSILON as BaseEPSILON;
use std::f64::{INFINITY, NEG_INFINITY};

use super::Shape;

use crate::intersection::IntersectionList;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

const EPSILON: f64 = BaseEPSILON * 1e6;

pub struct Cylinder {
    material: Material,
    inv_tf: Matrix,
    xpose_inv_tf: Matrix,
    min: f64,
    max: f64,
    closed: bool,
}

impl Cylinder {
    pub fn new(material: Material, tf: Matrix, bounds: (f64, f64), closed: bool) -> Self {
        let inv_tf = tf.inverse();

        Self {
            material,
            inv_tf: inv_tf.clone(),
            xpose_inv_tf: inv_tf.transpose(),
            min: bounds.0,
            max: bounds.1,
            closed,
        }
    }

    pub fn default() -> Self {
        Self {
            material: Material::default(),
            inv_tf: Matrix::eye(4),
            xpose_inv_tf: Matrix::eye(4),
            min: NEG_INFINITY,
            max: INFINITY,
            closed: false,
        }
    }

    pub fn with_material(mut self, m: Material) -> Self {
        self.material = m;
        self
    }

    pub fn with_transform(mut self, tf: Matrix) -> Self {
        let inv_tf = tf.inverse();

        self.inv_tf = inv_tf.clone();
        self.xpose_inv_tf = inv_tf.transpose();

        self
    }

    pub fn with_bounds(mut self, min: f64, max: f64, closed: bool) -> Self {
        self.min = min;
        self.max = max;
        self.closed = closed;

        self
    }

    fn intersect_caps(&self, r: Ray) -> Vec<f64> {
        let mut xs = Vec::new();

        // only check the caps intersection if the cylinder is
        // not closed and the direction has a y component.
        if !self.closed || r.direction.y.abs() < EPSILON {
            return xs;
        }

        // check intersection at the lower end cap
        let t = (self.min - r.origin.y) / r.direction.y;
        if check_cap(r, t) {
            xs.push(t);
        }

        // check intersection at he upper end cap
        let t = (self.max - r.origin.y) / r.direction.y;
        if check_cap(r, t) {
            xs.push(t);
        }

        xs
    }
}

fn check_cap(r: Ray, t: f64) -> bool {
    let x = r.origin.x + t * r.direction.x;
    let z = r.origin.z + t * r.direction.z;

    x * x + z * z <= 1.
}

impl Shape for Cylinder {
    #[inline(always)]
    fn inverse_transform(&self) -> &Matrix {
        &self.inv_tf
    }

    #[inline(always)]
    fn transpose_inverse(&self) -> &Matrix {
        &self.xpose_inv_tf
    }

    #[inline(always)]
    fn material(&self) -> &Material {
        &self.material
    }

    fn local_intersect(&self, r: Ray) -> Vec<f64> {
        let mut xs = self.intersect_caps(r);

        let a = r.direction.x * r.direction.x + r.direction.z * r.direction.z;
        if a.abs() < EPSILON {
            // if parallel to the y axis
            return xs;
        }

        let b = 2. * r.origin.x * r.direction.x + 2. * r.origin.z * r.direction.z;
        let c = r.origin.x * r.origin.x + r.origin.z * r.origin.z - 1.;

        let disc = b * b - 4. * a * c;
        if disc < 0. {
            // ray does not intersect cylinder
            return xs;
        }

        let mut t0 = (-b - disc.sqrt()) / (2. * a);
        let mut t1 = (-b + disc.sqrt()) / (2. * a);
        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }

        let y0 = r.origin.y + t0 * r.direction.y;
        if self.min < y0 && y0 < self.max {
            xs.push(t0);
        }

        let y1 = r.origin.y + t1 * r.direction.y;
        if self.min < y1 && y1 < self.max {
            xs.push(t1);
        }

        xs
    }

    fn local_normal_at(&self, p: Point) -> Vector {
        let dist = p.x * p.x + p.z * p.z;
        if dist < 1. && p.y >= self.max - EPSILON {
            return Vector::new(0., 1., 0.);
        } else if dist < 1. && p.y <= self.min + EPSILON {
            return Vector::new(0., -1., 0.);
        }

        Vector::new(p.x, 0., p.z)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn default() {
        let cyl = Cylinder::default();
        assert_eq!(cyl.min, NEG_INFINITY);
        assert_eq!(cyl.max, INFINITY);
        assert_eq!(cyl.closed, false);
    }

    mod local_intersect {
        use super::*;

        #[test]
        fn ray_misses_a_cylinder() {
            let cyl = Cylinder::default();

            let params = vec![
                (Point::new(1., 0., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0., 0., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0., 0., -5.), Vector::new(1., 1., 1.)),
            ];
            for (origin, direction) in params {
                let r = Ray::new(origin, direction.normalize());
                let xs = cyl.local_intersect(r);
                assert_eq!(xs.len(), 0);
            }
        }

        #[test]
        fn ray_hits_a_cylinder() {
            let cyl = Cylinder::default();

            let params = vec![
                (Point::new(1., 0., -5.), Vector::new(0., 0., 1.), 5., 5.),
                (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 4., 6.),
                (
                    Point::new(0.5, 0., -5.),
                    Vector::new(0.1, 1., 1.),
                    6.80798191702732,
                    7.088723439378861,
                ),
            ];
            for (origin, direction, t0, t1) in params {
                let r = Ray::new(origin, direction.normalize());
                let xs = cyl.local_intersect(r);
                assert_eq!(xs.len(), 2);
                assert_eq!(xs[0], t0);
                assert_eq!(xs[1], t1);
            }
        }

        #[test]
        fn intersecting_a_constrained_cylinder() {
            let cyl = Cylinder::default().with_bounds(1., 2., false);

            let params = vec![
                (Point::new(0., 1.5, 0.), Vector::new(0.1, 1., 0.), 0),
                (Point::new(0., 3., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 2., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 1., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 1.5, -2.), Vector::new(0., 0., 1.), 2),
            ];
            for (origin, direction, count) in params {
                let r = Ray::new(origin, direction.normalize());
                let xs = cyl.local_intersect(r);
                assert_eq!(xs.len(), count);
            }
        }

        #[test]
        fn intersecting_the_caps_of_a_closed_cylinder() {
            let cyl = Cylinder::default().with_bounds(1., 2., true);

            let params = vec![
                (Point::new(0., 3., 0.), Vector::new(0., -1., 0.), 2),
                (Point::new(0., 3., -2.), Vector::new(0., -1., 2.), 2),
                (Point::new(0., 4., -2.), Vector::new(0., -1., 1.), 2), // corner case
                (Point::new(0., 0., -2.), Vector::new(0., 1., 2.), 2),
                (Point::new(0., -1., -2.), Vector::new(0., 1., 1.), 2), // corner case
            ];

            for (origin, direction, count) in params {
                let r = Ray::new(origin, direction.normalize());
                let xs = cyl.local_intersect(r);
                assert_eq!(xs.len(), count);
            }
        }
    }

    mod local_normal_at {
        use super::*;

        #[test]
        fn cylinder_sides() {
            let cyl = Cylinder::default();

            let params = vec![
                (Point::new(1., 0., 0.), Vector::new(1., 0., 0.)),
                (Point::new(0., 5., -1.), Vector::new(0., 0., -1.)),
                (Point::new(0., -2., 1.), Vector::new(0., 0., 1.)),
                (Point::new(-1., 1., 0.), Vector::new(-1., 0., 0.)),
            ];
            for (point, normal) in params {
                let v = cyl.local_normal_at(point);
                assert_eq!(v, normal);
            }
        }

        #[test]
        fn cylinder_end_caps() {
            let cyl = Cylinder::default().with_bounds(1., 2., true);

            let params = vec![
                (Point::new(0., 1., 0.), Vector::new(0., -1., 0.)),
                (Point::new(0.5, 1., 0.), Vector::new(0., -1., 0.)),
                (Point::new(0., 1., 0.5), Vector::new(0., -1., 0.)),
                (Point::new(0., 2., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0.5, 2., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0., 2., 0.5), Vector::new(0., 1., 0.)),
            ];
            for (point, normal) in params {
                let v = cyl.local_normal_at(point);
                assert_eq!(v, normal);
            }
        }
    }
}
