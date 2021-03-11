use std::f64::INFINITY;
use std::mem;

use super::Shape;

use crate::intersection::IntersectionList;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

const EPSILON: f64 = 1e-9;

pub struct Cube {
    material: Material,
    inv_tf: Matrix,
    xpose_inv_tf: Matrix,
}

impl Cube {
    pub fn new() -> Self {
        Self {
            material: Material::default(),
            inv_tf: Matrix::eye(4),
            xpose_inv_tf: Matrix::eye(4),
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

    pub fn set_material(&mut self, m: Material) {
        self.material = m;
    }

    pub fn set_transform(&mut self, tf: Matrix) {
        let inv_tf = tf.inverse();

        self.inv_tf = inv_tf.clone();
        self.xpose_inv_tf = inv_tf.transpose()
    }
}

impl Shape for Cube {
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
        let (xtmin, xtmax) = check_axis(r.origin.x, r.direction.x);
        let (ytmin, ytmax) = check_axis(r.origin.y, r.direction.y);
        let (ztmin, ztmax) = check_axis(r.origin.z, r.direction.z);

        let tmin = f64::max(f64::max(xtmin, ytmin), ztmin);
        let tmax = f64::min(f64::min(xtmax, ytmax), ztmax);

        if tmin > tmax {
            return Vec::new(); // no intersections with ray and cube
        }

        vec![tmin, tmax]
    }
    fn local_normal_at(&self, p: Point) -> Vector {
        match f64::max(f64::max(p.x.abs(), p.y.abs()), p.z.abs()) {
            max if max == p.x.abs() => Vector::new(p.x, 0., 0.),
            max if max == p.y.abs() => Vector::new(0., p.y, 0.),
            _ => Vector::new(0., 0., p.z), // max == p.z.abs()
        }
    }
}

fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
    let tmin_numerator = -1. - origin;
    let tmax_numerator = 1. - origin;

    let mut tmin: f64;
    let mut tmax: f64;
    if direction.abs() > EPSILON {
        tmin = tmin_numerator / direction;
        tmax = tmax_numerator / direction;
    } else {
        tmin = tmin_numerator * INFINITY;
        tmax = tmax_numerator * INFINITY;
    }

    if tmin > tmax {
        mem::swap(&mut tmin, &mut tmax);
    }

    (tmin, tmax)
}

#[cfg(test)]
mod test {
    use super::*;

    mod local_intersect {
        use super::*;

        #[test]
        fn ray_intersects_a_cube() {
            let c = Cube::new();
            let params = vec![
                (Point::new(5., 0.5, 0.), Vector::new(-1., 0., 0.), 4., 6.), // +x
                (Point::new(-5., 0.5, 0.), Vector::new(1., 0., 0.), 4., 6.), // -x
                (Point::new(0.5, 5., 0.), Vector::new(0., -1., 0.), 4., 6.), // +y
                (Point::new(0.5, -5., 0.), Vector::new(0., 1., 0.), 4., 6.), // -y
                (Point::new(0.5, 0., 5.), Vector::new(0., 0., -1.), 4., 6.), // +z
                (Point::new(0.5, 0., -5.), Vector::new(0., 0., 1.), 4., 6.), // -z
            ];

            for (origin, direction, t1, t2) in params {
                let r = Ray::new(origin, direction);
                let xs = c.local_intersect(r);
                assert_eq!(xs.len(), 2);
                assert_eq!(xs[0], t1);
                assert_eq!(xs[1], t2);
            }
        }

        #[test]
        fn ray_misses_a_cube() {
            let c = Cube::new();
            let params = vec![
                (Point::new(-2., 0., 0.), Vector::new(0.2673, 0.5345, 0.8018)),
                (Point::new(0., -2., 0.), Vector::new(0.8018, 0.2673, 0.5345)),
                (Point::new(0., 0., -2.), Vector::new(0.5345, 0.8018, 0.2673)),
                (Point::new(2., 0., 2.), Vector::new(0., 0., -1.)),
                (Point::new(0., 2., 2.), Vector::new(0., -1., 0.)),
                (Point::new(2., 2., 0.), Vector::new(-1., 0., 0.)),
            ];

            for (origin, direction) in params {
                let r = Ray::new(origin, direction);
                let xs = c.local_intersect(r);
                assert_eq!(xs.len(), 0);
            }
        }
    }

    #[test]
    fn local_normal_at() {
        let c = Cube::new();

        let params = vec![
            (Point::new(1., 0.5, -0.8), Vector::new(1., 0., 0.)),
            (Point::new(-1., -0.2, 0.9), Vector::new(-1., 0., 0.)),
            (Point::new(-0.4, 1., -0.1), Vector::new(0., 1., 0.)),
            (Point::new(0.3, -1., -0.7), Vector::new(0., -1., 0.)),
            (Point::new(-0.6, 0.3, 1.), Vector::new(0., 0., 1.)),
            (Point::new(0.4, 0.4, -1.), Vector::new(0., 0., -1.)),
            // the normal at corners will be in the x-direction
            (Point::new(1., 1., 1.), Vector::new(1., 0., 0.)),
            (Point::new(-1., -1., -1.), Vector::new(-1., 0., 0.)),
        ];

        for (point, normal) in params {
            assert_eq!(c.local_normal_at(point), normal)
        }
    }
}
