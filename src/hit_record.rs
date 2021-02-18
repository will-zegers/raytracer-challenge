use crate::intersection::Intersection;
use crate::point::Point;
use crate::ray::Ray;
use crate::sphere::Sphere;
use crate::vector::Vector;

pub struct HitRecord<'a> {
    pub t: f64,
    pub object: &'a Sphere,
    pub point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub inside: bool,
}

#[allow(dead_code)]
impl<'a> HitRecord<'a> {
    pub fn new(i: &Intersection<'a>, r: &Ray) -> Self {
        let t = i.t;
        let object = i.object;
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
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new() {
        // precomputing the state of an intersection, when an intersection occurs on the outside
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let i = Intersection::new(4., &s);
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
        let i = Intersection::new(1., &s);
        let h = HitRecord::new(&i, &r);

        assert_eq!(h.point, Point::new(0., 0., 1.));
        assert_eq!(h.eyev, Vector::new(0., 0., -1.));
        assert!(h.inside);
        assert_eq!(h.normalv, Vector::new(0., 0., -1.))
    }
}
