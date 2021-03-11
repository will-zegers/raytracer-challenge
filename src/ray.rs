use std::rc::Rc;

use crate::intersection::{Intersection, IntersectionList};
use crate::matrix::Matrix;
use crate::point::Point;
use crate::shape::Shape;
use crate::vector::Vector;

#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vector,
}

impl Ray {
    pub fn new(origin: Point, direction: Vector) -> Self {
        Self { origin, direction }
    }

    pub fn at(&self, t: f64) -> Point {
        self.origin + t * self.direction
    }

    pub fn intersects(&self, s: Rc<dyn Shape>) -> IntersectionList {
        let ray = Self::transform(self, &s.inverse_transform());

        let mut xs = IntersectionList::new();
        for hit in s.local_intersect(ray) {
            xs.push(Intersection::new(hit, s.clone()));
        }

        xs
    }

    pub fn transform(ray: &Self, t: &Matrix) -> Self {
        Self {
            origin: t * ray.origin,
            direction: t * ray.direction,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::shape::Sphere;

    #[test]
    fn new() {
        let origin = Point::new(1., 2., 3.);
        let direction = Vector::new(4., 5., 6.);
        let r = Ray::new(origin, direction);

        assert_eq!(r.origin, origin);
        assert_eq!(r.direction, direction);
    }

    #[test]
    fn at() {
        let r = Ray::new(Point::new(2., 3., 4.), Vector::new(1., 0., 0.));
        assert_eq!(r.at(0.), Point::new(2., 3., 4.));
        assert_eq!(r.at(1.), Point::new(3., 3., 4.));
        assert_eq!(r.at(-1.), Point::new(1., 3., 4.));
        assert_eq!(r.at(2.5), Point::new(4.5, 3., 4.));
    }

    #[test]
    fn intersects() {
        // ray goes through the unit sphere's center
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = r.intersects(Rc::new(s));
        assert_eq!(xs[0].t, 4.);
        assert_eq!(xs[1].t, 6.);

        // ray is tangential to unit sphere
        let r = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = r.intersects(Rc::new(s));
        assert_eq!(xs[0].t, 5.);
        assert_eq!(xs[1].t, 5.);

        // ray misses a sphere
        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = r.intersects(Rc::new(s));
        assert!(xs.is_empty());

        // ray's origin is inside the sphere
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = r.intersects(Rc::new(s));
        assert_eq!(xs[0].t, -1.);
        assert_eq!(xs[1].t, 1.);

        // sphere is behind the ray
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = r.intersects(Rc::new(s));
        assert_eq!(xs[0].t, -6.);
        assert_eq!(xs[1].t, -4.);
    }

    #[test]
    fn transform() {
        // translation
        let r = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        let m = Matrix::translation(3., 4., 5.);
        let r2 = Ray::transform(&r, &m);
        assert_eq!(r2.origin, Point::new(4., 6., 8.));
        assert_eq!(r2.direction, Vector::new(0., 1., 0.));

        // scaling
        let r = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        let m = Matrix::scaling(2., 3., 4.);
        let r2 = Ray::transform(&r, &m);
        assert_eq!(r2.origin, Point::new(2., 6., 12.));
        assert_eq!(r2.direction, Vector::new(0., 3., 0.));

        // intersecting a scaled sphere with a ray
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new().set_transform(Matrix::scaling(2., 2., 2.));
        let xs = r.intersects(Rc::new(s));
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 3.);
        assert_eq!(xs[1].t, 7.);

        // intersecting a translated sphere with a ray
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new().set_transform(Matrix::translation(5., 0., 0.));
        assert!(r.intersects(Rc::new(s)).is_empty());
    }
}
