use std::f64::{EPSILON as BaseEPSILON, INFINITY, NEG_INFINITY};

use super::Shape;

use crate::intersection::IntersectionList;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

const EPSILON: f64 = BaseEPSILON * 1e6;

pub struct Cone {
    material: Material,
    inv_tf: Matrix,
    xpose_inv_tf: Matrix,
    min: f64,
    max: f64,
    closed: bool,
}

impl Cone {
    pub fn new(material: Material, tf: Matrix, bounds: (f64, f64), closed: bool) -> Self {
        let inv_tf = tf.inverse();
        Self {
            material,
            inv_tf: inv_tf.inverse().clone(),
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

    fn intersects_caps(&self, r: Ray) -> Vec<f64> {
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
    let y = r.origin.y + t * r.direction.y;
    let z = r.origin.z + t * r.direction.z;

    x * x + z * z <= y * y
}

impl Shape for Cone {
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
        let mut xs = self.intersects_caps(r);

        let a = r.direction.x * r.direction.x - r.direction.y * r.direction.y
            + r.direction.z * r.direction.z;
        let b = 2. * r.origin.x * r.direction.x - 2. * r.origin.y * r.direction.y
            + 2. * r.origin.z * r.direction.z;

        // no intersection
        if a.abs() < EPSILON && b.abs() < EPSILON {
            return xs;
        }

        let c = r.origin.x * r.origin.x - r.origin.y * r.origin.y + r.origin.z * r.origin.z;

        // ray is parallel to one half, but intersects the other
        if a.abs() < EPSILON {
            xs.push(-c / (2. * b));
            return xs;
        }

        let disc = b * b - 4. * a * c;
        if disc < 0. {
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
        if p.y >= self.max - EPSILON {
            return Vector::new(0., 1., 0.);
        } else if p.y <= self.min + EPSILON {
            return Vector::new(0., -1., 0.);
        }

        let y = if p.y > 0. {
            -(p.x * p.x + p.z * p.z).sqrt()
        } else {
            (p.x * p.x + p.z * p.z).sqrt()
        };
        Vector::new(p.x, y, p.z)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    mod local_intersect {
        use super::*;

        #[test]
        fn intersecting_a_cone_with_a_ray() {
            let shape = Cone::default();

            let params = vec![
                (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 5., 5.),
                (
                    Point::new(0., 0., -5.),
                    Vector::new(1., 1., 1.),
                    8.660254037844386,
                    8.660254037844386,
                ),
                (
                    Point::new(1., 1., -5.),
                    Vector::new(-0.5, -1., 1.),
                    4.550055679356349,
                    49.449944320643645,
                ),
            ];
            for (origin, direction, t0, t1) in params {
                let r = Ray::new(origin, direction.normalize());
                let xs = shape.local_intersect(r);
                assert_eq!(xs.len(), 2);
                assert_eq!(xs[0], t0);
                assert_eq!(xs[1], t1);
            }
        }

        #[test]
        fn intersecting_a_cone_with_a_ray_parallel_to_one_of_its_halves() {
            let shape = Cone::default();

            let direction = Vector::new(0., 1., 1.).normalize();
            let r = Ray::new(Point::new(0., 0., -1.), direction);
            let xs = shape.local_intersect(r);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0], 0.3535533905932738);
        }

        #[test]
        fn intersection_cones_end_caps() {
            let shape = Cone::default()
                .with_bounds(-0.5, 0.5, true);

            let params = vec![
                (Point::new(0., 0., -5.), Vector::new(0., 1., 0.), 0),
                (Point::new(0., 0., -0.25), Vector::new(0., 1., 1.), 2),
                (Point::new(0., 0., 0.25), Vector::new(0., 1., 0.), 4),
            ];
            for (origin, direction, count) in params {
                let r = Ray::new(origin, direction.normalize());
                let xs = shape.local_intersect(r);
                assert_eq!(xs.len(), count);
            }
        }
    }

    #[test]
    fn local_normal_at() {
        let shape = Cone::default();

        let params = vec![
            (Point::new(0., 0., 0.), Vector::new(0., 0., 0.)),
            (Point::new(1., 1., 1.), Vector::new(1., -2_f64.sqrt(), 1.)),
            (Point::new(-1., -1., 0.), Vector::new(-1., 1., 0.)),
        ];
        for (point, normal) in params {
            assert_eq!(shape.local_normal_at(point), normal);
        }
    }
}
