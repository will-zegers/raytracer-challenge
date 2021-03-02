use std::ops::Index;
use std::rc::Rc;

use crate::shape::{Shape, Sphere};

#[cfg_attr(test, derive(Clone, Debug))]
pub struct Intersection {
    pub t: f64,
    pub object: Rc<dyn Shape>,
}

impl Intersection {
    pub fn new(t: f64, object: Rc<dyn Shape>) -> Self {
        Self { t, object }
    }
}

#[cfg_attr(test, derive(Debug))]
pub struct IntersectionList {
    objects: Vec<Intersection>,
}

impl IntersectionList {
    pub fn empty() -> Self {
        Self {
            objects: Vec::new(),
        }
    }

    #[cfg(test)]
    pub fn new(mut intersections: Vec<Intersection>) -> Self {
        IntersectionList::sort(&mut intersections);
        Self {
            objects: intersections,
        }
    }

    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.objects.len()
    }

    pub fn is_empty(&self) -> bool {
        self.objects.is_empty()
    }

    pub fn hit(&self) -> Option<&Intersection> {
        self.objects.iter().find(|&x| x.t >= 0.)
    }

    pub fn push(&mut self, x: Intersection) {
        self.objects.push(x);
    }

    pub fn extend(&mut self, xs: Self) {
        self.objects.extend(xs.objects);
        IntersectionList::sort(&mut self.objects);
    }

    fn sort(xs: &mut Vec<Intersection>) {
        xs.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    }
}

impl Index<usize> for IntersectionList {
    type Output = Intersection;

    fn index(&self, i: usize) -> &Self::Output {
        &self.objects[i]
    }
}

#[cfg(test)]
mod test {
    use super::*;

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

        let xs = IntersectionList::new(vec![i1, i2]);

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
        let xs = IntersectionList::new(vec![i2, i1.clone()]);
        let i = xs.hit().unwrap();
        assert_eq!(*i, i1);

        // the hit, when some intersections have a negative t
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(-1., s.clone());
        let i2 = Intersection::new(1., s.clone());
        let xs = IntersectionList::new(vec![i2.clone(), i1]);
        let i = xs.hit().unwrap();
        assert_eq!(*i, i2);

        // the hit, when all intersections have a negative t
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(-2., s.clone());
        let i2 = Intersection::new(-1., s.clone());
        let xs = IntersectionList::new(vec![i2, i1]);
        let i = xs.hit();
        assert!(i.is_none());

        // the hit, always the lowest non-negative intersection
        let s = Rc::new(Sphere::new());
        let i1 = Intersection::new(5., s.clone());
        let i2 = Intersection::new(7., s.clone());
        let i3 = Intersection::new(-3., s.clone());
        let i4 = Intersection::new(2., s.clone());
        let xs = IntersectionList::new(vec![i1, i2, i3, i4.clone()]);
        let i = xs.hit().unwrap();
        assert_eq!(*i, i4);
    }
}
