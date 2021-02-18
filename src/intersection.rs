use std::ops::Index;

use crate::sphere::Sphere;

#[derive(Clone, Debug)]
pub struct Intersection<'a> {
    pub t: f64,
    pub object: &'a Sphere,
}

impl<'a> Intersection<'a> {
    pub fn new(t: f64, object: &'a Sphere) -> Self {
        Self { t, object }
    }
}

#[derive(Debug)]
pub struct IntersectionList<'a> {
    objects: Vec<Intersection<'a>>,
}

impl<'a> IntersectionList<'a> {
    pub fn empty() -> Self {
        Self {
            objects: Vec::new(),
        }
    }

    pub fn new(mut intersections: Vec<Intersection<'a>>) -> Self {
        IntersectionList::sort(&mut intersections);
        Self {
            objects: intersections,
        }
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.objects.len()
    }

    #[allow(dead_code)]
    pub fn hit(&self) -> Option<&Intersection> {
        self.objects.iter().find(|&x| x.t >= 0.)
    }

    #[allow(dead_code)]
    pub fn extend(&mut self, xs: Self) {
        self.objects.extend(xs.objects);
        IntersectionList::sort(&mut self.objects);
    }

    fn sort(xs: &mut Vec<Intersection<'a>>) {
        xs.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    }
}

impl<'a> Index<usize> for IntersectionList<'a> {
    type Output = Intersection<'a>;

    fn index(&self, i: usize) -> &Self::Output {
        &self.objects[i]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    impl<'a> PartialEq<Intersection<'a>> for Intersection<'a> {
        fn eq(&self, rhs: &Self) -> bool {
            self.t == rhs.t && self.object == rhs.object
        }
    }

    #[test]
    fn new() {
        let s = Sphere::default();
        let i = Intersection::new(3.5, &s);
        assert_eq!(i.t, 3.5);
        assert_eq!(i.object, &s);

        let s = Sphere::default();
        let i1 = Intersection::new(1., &s);
        let i2 = Intersection::new(2., &s);

        let xs = IntersectionList::new(vec![i1, i2]);

        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit() {
        // the hit, when all intersections have a positive t
        let s = Sphere::default();
        let i1 = Intersection::new(1., &s);
        let i2 = Intersection::new(2., &s);
        let xs = IntersectionList::new(vec![i2, i1.clone()]);
        let i = xs.hit().unwrap();
        assert_eq!(*i, i1);

        // the hit, when some intersections have a negative t
        let s = Sphere::default();
        let i1 = Intersection::new(-1., &s);
        let i2 = Intersection::new(1., &s);
        let xs = IntersectionList::new(vec![i2.clone(), i1]);
        let i = xs.hit().unwrap();
        assert_eq!(*i, i2);

        // the hit, when all intersections have a negative t
        let s = Sphere::default();
        let i1 = Intersection::new(-2., &s);
        let i2 = Intersection::new(-1., &s);
        let xs = IntersectionList::new(vec![i2, i1]);
        let i = xs.hit();
        assert!(i.is_none());

        // the hit, always the lowest non-negative intersection
        let s = Sphere::default();
        let i1 = Intersection::new(5., &s);
        let i2 = Intersection::new(7., &s);
        let i3 = Intersection::new(-3., &s);
        let i4 = Intersection::new(2., &s);
        let xs = IntersectionList::new(vec![i1, i2, i3, i4.clone()]);
        let i = xs.hit().unwrap();
        assert_eq!(*i, i4);
    }
}
