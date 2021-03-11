use std::f64::EPSILON as BaseEPSILON;
use std::ptr;
use std::rc::Rc;

use crate::intersection::{Intersection, IntersectionList};
use crate::point::Point;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::vector::Vector;

const EPSILON: f64 = BaseEPSILON * 1e6;

pub struct HitRecord {
    pub t: f64,
    pub object: Rc<dyn Shape>,
    pub point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub inside: bool,
    pub over_point: Point,
    pub under_point: Point,
    pub reflectv: Vector,
    pub n1: f64,
    pub n2: f64,
}

impl HitRecord {
    pub fn new(i: &Intersection, r: &Ray, xs: &IntersectionList) -> Self {
        let mut containers: Vec<Rc<dyn Shape>> = Vec::new();
        let mut n1 = 1.;
        let mut n2 = 1.;

        'refraction_indices: for j in 0..xs.len() {
            if ptr::eq(&i.object, &xs[j].object) {
                if let Some(obj) = containers.last() {
                    n1 = obj.material().refractive_index;
                }
            }
            match containers
                .iter()
                .position(|obj| Rc::ptr_eq(&obj, &xs[j].object))
            {
                Some(found) => {
                    containers.remove(found);
                }
                None => containers.push(xs[j].object.clone()),
            }

            if ptr::eq(&i.object, &xs[j].object) {
                if let Some(obj) = containers.last() {
                    n2 = obj.material().refractive_index;
                }
                break 'refraction_indices;
            }
        }

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
            under_point: point - EPSILON * normalv,
            reflectv: r.direction.reflect(&normalv),
            n1,
            n2,
        }
    }

    pub fn schlick(&self) -> f64 {
        // find the cosine of the angle between the eye and normal vectors
        let mut cos = Vector::dot(&self.eyev, &self.normalv);

        if self.n1 > self.n2 {
            let n = self.n1 / self.n2;
            let sin2_t = n * n * (1. - cos * cos);
            if sin2_t > 1. {
                return 1.; // total internal reflection
            }

            let cos_t = (1. - sin2_t).sqrt();
            cos = cos_t;
        }

        let mut r0 = (self.n1 - self.n2) / (self.n1 + self.n2);
        r0 *= r0; // r0^2

        r0 + (1. - r0) * (1. - cos).powi(5)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::light::PointLight;
    use crate::material::Material;
    use crate::matrix::Matrix;
    use crate::shape::{Plane, Sphere};
    use crate::world::World;

    #[test]
    fn new() {
        // precomputing the state of an intersection, when an intersection occurs on the outside
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let i = Intersection::new(4., Rc::new(s));
        let h = HitRecord::new(&i, &r, &IntersectionList::new());

        assert_eq!(h.t, i.t);
        // assert_eq!(h.object, i.object);
        assert_eq!(h.point, Point::new(0., 0., -1.));
        assert_eq!(h.eyev, Vector::new(0., 0., -1.));
        assert_eq!(h.normalv, Vector::new(0., 0., -1.));
        assert!(!h.inside);

        // the hit, when the intersection occurs on the inside
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let i = Intersection::new(1., Rc::new(s));
        let h = HitRecord::new(&i, &r, &IntersectionList::new());

        assert_eq!(h.point, Point::new(0., 0., 1.));
        assert_eq!(h.eyev, Vector::new(0., 0., -1.));
        assert!(h.inside);
        assert_eq!(h.normalv, Vector::new(0., 0., -1.))
    }

    #[test]
    fn over_point() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let t = Matrix::translation(0., 0., 1.);
        let shape = Sphere::new().set_transform(t);
        let i = Intersection::new(5., Rc::new(shape));

        let rec = HitRecord::new(&i, &r, &IntersectionList::new());
        assert!(rec.over_point.z < EPSILON / 2.);
        assert!(rec.point.z > rec.over_point.z);
    }

    #[test]
    fn reflectv() {
        // precomputing the reflection vector, reflectv
        let shape = Rc::new(Plane::default());
        let r = Ray::new(
            Point::new(0., 1., -1.),
            Vector::new(0., -f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.),
        );
        let xs = r.intersects(shape.clone());
        let i = Intersection::new(xs[0].t, shape);
        let rec = HitRecord::new(&i, &r, &IntersectionList::new());
        assert_eq!(
            rec.reflectv,
            Vector::new(0., f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.)
        );
    }

    #[test]
    fn refraction() {
        let mut m_a = Material::default();
        m_a.refractive_index = 1.5;
        let a = Sphere::glass()
            .set_transform(Matrix::scaling(2., 2., 2.))
            .set_material(m_a);
        let a = Rc::new(a);

        let mut m_b = Material::default();
        m_b.refractive_index = 2.0;
        let b = Sphere::glass()
            .set_transform(Matrix::translation(0., 0., -0.25))
            .set_material(m_b);
        let b = Rc::new(b);

        let mut m_c = Material::default();
        m_c.refractive_index = 2.5;
        let c = Sphere::glass()
            .set_transform(Matrix::translation(0., 0., 0.25))
            .set_material(m_c);
        let c = Rc::new(c);

        let r = Ray::new(Point::new(0., 0., -4.), Vector::new(0., 0., 1.));
        let xs = vec![
            Intersection::new(2., a.clone()),
            Intersection::new(2.75, b.clone()),
            Intersection::new(3.25, c.clone()),
            Intersection::new(4.25, b.clone()),
            Intersection::new(5.25, c.clone()),
            Intersection::new(6., a.clone()),
        ];

        let rec = HitRecord::new(&xs[0], &r, &xs);
        assert_eq!(rec.n1, 1.0);
        assert_eq!(rec.n2, 1.5);

        let rec = HitRecord::new(&xs[1], &r, &xs);
        assert_eq!(rec.n1, 1.5);
        assert_eq!(rec.n2, 2.0);

        let rec = HitRecord::new(&xs[2], &r, &xs);
        assert_eq!(rec.n1, 2.0);
        assert_eq!(rec.n2, 2.5);

        let rec = HitRecord::new(&xs[3], &r, &xs);
        assert_eq!(rec.n1, 2.5);
        assert_eq!(rec.n2, 2.5);

        let rec = HitRecord::new(&xs[4], &r, &xs);
        assert_eq!(rec.n1, 2.5);
        assert_eq!(rec.n2, 1.5);

        let rec = HitRecord::new(&xs[5], &r, &xs);
        assert_eq!(rec.n1, 1.5);
        assert_eq!(rec.n2, 1.0);
    }

    #[test]
    fn under_point() {
        // the under point is the offset below the surface
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::glass().set_transform(Matrix::translation(0., 0., 1.));
        let s = Rc::new(s);
        let i = Intersection::new(5., s);
        let xs = vec![i.clone()];

        let rec = HitRecord::new(&i, &r, &xs);

        assert!(rec.under_point.z > EPSILON / 2.);
        assert!(rec.point.z < rec.under_point.z);
    }

    mod schlick {
        use super::*;

        use crate::light::PointLight;
        use crate::world::World;

        const TOL: f64 = 1e-9;

        #[test]
        fn total_internal_refraction() {
            let shape = Rc::new(Sphere::glass());
            let r = Ray::new(
                Point::new(0., 0., f64::sqrt(2.) / 2.),
                Vector::new(0., 1., 0.),
            );
            let mut xs: Vec<Intersection> = Vec::new();
            for t in shape.local_intersect(r) {
                xs.push(Intersection::new(t, shape.clone()));
            }
            let rec = HitRecord::new(&xs[1], &r, &xs);
            assert_eq!(rec.schlick(), 1.);
        }

        #[test]
        fn perpendicular_viewing_angle() {
            let shape = Rc::new(Sphere::glass());
            let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
            let mut xs: Vec<Intersection> = Vec::new();
            for t in shape.local_intersect(r) {
                xs.push(Intersection::new(t, shape.clone()));
            }
            let rec = HitRecord::new(&xs[1], &r, &xs);
            assert!((rec.schlick() - 0.04).abs() < TOL);
        }

        #[test]
        fn small_angle_and_n2_gt_n1() {
            let shape = Rc::new(Sphere::glass());
            let r = Ray::new(Point::new(0., 0.99, -2.), Vector::new(0., 0., 1.));
            let mut xs: Vec<Intersection> = Vec::new();
            for t in shape.local_intersect(r) {
                xs.push(Intersection::new(t, shape.clone()));
            }
            let rec = HitRecord::new(&xs[0], &r, &xs);
            assert!((rec.schlick() - 0.488814383).abs() < TOL);
        }
    }
}
