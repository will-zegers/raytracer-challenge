use std::ops::Index;
use std::rc::Rc;

use crate::shape::{Shape, Sphere};

#[derive(Clone)]
#[cfg_attr(test, derive(Debug))]
pub struct Intersection {
    pub t: f64,
    pub object: Rc<dyn Shape>,
}

impl Intersection {
    pub fn new(t: f64, object: Rc<dyn Shape>) -> Self {
        Self { t, object }
    }
}

pub fn hit(xs: &IntersectionList) -> Option<&Intersection> {
    xs.iter().find(|&x| x.t >= 0.)
}

pub fn sort(xs: &mut IntersectionList) {
    xs.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
}

pub type IntersectionList = Vec<Intersection>;

#[cfg(test)]
mod test {
    use super::*;

    use crate::intersection;

    impl PartialEq<Intersection> for Intersection {
        fn eq(&self, rhs: &Self) -> bool {
            self.t == rhs.t // && self.object == rhs.object
        }
    }

    #[test]
    fn new() {
        let s = Rc::new(Sphere::new());
        let i = Intersection::new(3.5, s.clone());
        assert_eq!(i.t, 3.5);
        // assert_eq!(i.object, s.clone());

        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(1., s.clone());
        let i2 = Intersection::new(2., s.clone());

        let xs = vec![i1, i2];

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit() {
        // the hit, when all intersections have a positive t
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(1., s.clone());
        let i2 = Intersection::new(2., s.clone());
        let mut xs = vec![i2, i1.clone()];
        intersection::sort(&mut xs);
        let i = intersection::hit(&mut xs).unwrap();
        assert_eq!(*i, i1);

        // the hit, when some intersections have a negative t
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(-1., s.clone());
        let i2 = Intersection::new(1., s.clone());
        let mut xs = vec![i2.clone(), i1];
        intersection::sort(&mut xs);
        let i = intersection::hit(&mut xs).unwrap();
        assert_eq!(*i, i2);

        // the hit, when all intersections have a negative t
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(-2., s.clone());
        let i2 = Intersection::new(-1., s.clone());
        let mut xs = vec![i2, i1];
        intersection::sort(&mut xs);
        let i = intersection::hit(&mut xs);
        assert!(i.is_none());

        // the hit, always the lowest non-negative intersection
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(5., s.clone());
        let i2 = Intersection::new(7., s.clone());
        let i3 = Intersection::new(-3., s.clone());
        let i4 = Intersection::new(2., s.clone());
        let mut xs = vec![i1, i2, i3, i4.clone()];
        intersection::sort(&mut xs);
        let i = intersection::hit(&mut xs).unwrap();
        assert_eq!(*i, i4);
    }
}
