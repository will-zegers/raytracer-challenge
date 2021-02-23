use std::f64::EPSILON as BaseEPSILON;
use std::rc::Rc;

use crate::intersection::Intersection;
use crate::point::Point;
use crate::ray::Ray;
use crate::shape::Sphere;
use crate::vector::Vector;

const EPSILON: f64 = BaseEPSILON * 1e6;

pub struct HitRecord {
    pub t: f64,
    pub object: Rc<Sphere>,
    pub point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub inside: bool,
    pub over_point: Point,
}

#[allow(dead_code)]
impl HitRecord {
    pub fn new(i: &Intersection, r: &Ray) -> Self {
        let t = i.t;
        let object = i.object.clone();
        let point = r.at(t);
        let eyev = -r.direction;
        let mut normalv = object.normal_at(point);

        let inside;
        if Vector::dot(&normalv, &eyev) < 0. {
            inside = true;
            normalv = -normalv;
        } else {
            inside = false;
        }

        Self {
            t,
            object,
            point,
            eyev,
            normalv,
            inside,
            over_point: point + EPSILON * normalv,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::material::Material;
    use crate::matrix::Matrix;

    #[test]
    fn new() {
        // precomputing the state of an intersection, when an intersection occurs on the outside
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let i = Intersection::new(4., Rc::new(s));
        let h = HitRecord::new(&i, &r);

        assert_eq!(h.t, i.t);
        assert_eq!(h.object, i.object);
        assert_eq!(h.point, Point::new(0., 0., -1.));
        assert_eq!(h.eyev, Vector::new(0., 0., -1.));
        assert_eq!(h.normalv, Vector::new(0., 0., -1.));
        assert!(!h.inside);

        // the hit, when the intersection occurs on the inside
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let i = Intersection::new(1., Rc::new(s));
        let h = HitRecord::new(&i, &r);

        assert_eq!(h.point, Point::new(0., 0., 1.));
        assert_eq!(h.eyev, Vector::new(0., 0., -1.));
        assert!(h.inside);
        assert_eq!(h.normalv, Vector::new(0., 0., -1.))
    }

    #[test]
    fn over_point() {
        let r = Ray::new(Point::new(0., 0., -5.) , Vector::new(0., 0., 1.));
        let t = Matrix::translation(0., 0., 1.);
        let shape = Sphere::new(t, Material::default());
        let i = Intersection::new(5., Rc::new(shape));

        let rec = HitRecord::new(&i, &r);
        assert!(rec.over_point.z < EPSILON / 2.);
        assert!(rec.point.z > rec.over_point.z);
    }
}
