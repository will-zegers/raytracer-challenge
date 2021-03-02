use std::rc::Rc;

use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::color::Color;
use crate::hit_record::HitRecord;
use crate::intersection::IntersectionList;
use crate::light;
use crate::light::PointLight;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::patterns::Solid;
use crate::point::Point;
use crate::ray::Ray;
use crate::shape::{Shape, Sphere};

pub struct World {
    objects: Vec<Rc<dyn Shape>>,
    pub light: PointLight,
}

impl World {
    pub fn new(objects: Vec<Rc<dyn Shape>>, light: PointLight) -> Self {
        Self { objects, light }
    }

    #[cfg(test)]
    pub fn default() -> Self {
        let light = PointLight::default();
        let m1 = Material::new(Solid::new(Color::new(0.8, 1.0, 0.6)), 0.1, 0.7, 0.2, 200.);
        let s1 = Rc::new(Sphere::new().set_transform(Matrix::eye(4)).set_material(m1));
        let s2 = Rc::new(
            Sphere::new()
                .set_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .set_material(Material::default()),
        );

        Self {
            objects: vec![s1, s2],
            light,
        }
    }

    fn intersects(&self, r: &Ray) -> IntersectionList {
        let mut xs = IntersectionList::empty();
        for obj in &self.objects {
            let hits = r.intersects(obj.clone());
            if !hits.is_empty() {
                xs.extend(hits);
            }
        }
        xs
    }

    fn shade(&self, rec: &HitRecord) -> Color {
        light::lighting(
            rec.object.clone(),
            &self.light,
            &rec.over_point,
            &rec.eyev,
            &rec.normalv,
            self.is_shadowed(&rec.over_point),
        )
    }

    fn color_at(&self, r: &Ray) -> Color {
        let xs = self.intersects(r);
        return match xs.hit() {
            Some(x) => {
                let rec = HitRecord::new(xs.hit().unwrap(), r);
                self.shade(&rec)
            }
            None => Color::new(0., 0., 0.),
        };
    }

    pub fn render(&self, c: Camera) -> Canvas {
        let mut image = Canvas::new(c.hsize, c.vsize);
        for y in 0..c.vsize {
            for x in 0..c.hsize {
                let ray = c.pixel_to_ray(x, y);
                let color = self.color_at(&ray);
                image.write_pixel(x, y, color);
            }
        }
        image
    }

    fn is_shadowed(&self, p: &Point) -> bool {
        let v = self.light.position - *p;
        let distance = v.length();
        let direction = v.normalize();

        let r = Ray::new(*p, direction);
        let xs = self.intersects(&r);

        return match xs.hit() {
            Some(hit) => hit.t < distance,
            None => false,
        };
    }
}

#[cfg(test)]
mod test {
    use std::rc::Rc;

    use super::*;

    use crate::intersection::Intersection;
    use crate::point::Point;
    use crate::vector::Vector;

    #[test]
    fn default() {
        let m1 = Material::new(Solid::new(Color::new(0.8, 1.0, 0.6)), 0.1, 0.7, 0.2, 200.);
        let s1 = Rc::new(Sphere::new().set_transform(Matrix::eye(4)).set_material(m1));
        let s2 = Rc::new(
            Sphere::new()
                .set_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .set_material(Material::default()),
        );

        let w = World::default();
        assert_eq!(w.light, PointLight::default());
        // assert!(w.objects.iter().any(|s| **s == *s1));
        // assert!(w.objects.iter().any(|s| **s == *s2));
    }

    #[test]
    fn intersects() {
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = w.intersects(&r);
        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].t, 4.);
        assert_eq!(xs[1].t, 4.5);
        assert_eq!(xs[2].t, 5.5);
        assert_eq!(xs[3].t, 6.);
    }

    #[test]
    fn shade() {
        // shading an intersection
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = w.intersects(&r);
        assert_ne!(xs.len(), 0);
        assert_eq!(xs[0], Intersection::new(4.0, w.objects[0].clone()));

        let rec = HitRecord::new(&xs[0], &r);
        let c = w.shade(&rec);
        assert_eq!(c, Color::new(0.380661193, 0.475826491, 0.285495894));

        // shading an intersection from the inside
        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0.25, 0.), Color::new(1., 1., 1.));
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let xs = w.intersects(&r);
        let rec = HitRecord::new(&xs[2], &r);
        let c = w.shade(&rec);
        assert_eq!(c, Color::new(0.904984472, 0.904984472, 0.904984472));

        // shade() is given an intersection in shadow
        let s1 = Rc::new(Sphere::new());
        let t2 = Matrix::translation(0., 0., 10.);
        let s2 = Rc::new(
            Sphere::new()
                .set_transform(t2)
                .set_material(Material::default()),
        );
        let light = PointLight::new(Point::new(0., 0., -10.), Color::new(1., 1., 1.));
        let w = World::new(vec![s1, s2], light);

        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let xs = w.intersects(&r);
        let hit = xs.hit().unwrap();
        let rec = HitRecord::new(&hit, &r);
        let c = w.shade(&rec);
        assert_eq!(c, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn color_at() {
        // the color when a ray misses
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0., 0., 0.));

        // the color when a ray hits
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0.380661193, 0.475826491, 0.285495894));

        // the color with an intersection behind the ray
        let light = PointLight::default();

        let m1 = Material::new(Solid::new(Color::new(0.8, 1.0, 0.6)), 1., 0.7, 0.2, 200.);
        let s1 = Rc::new(Sphere::new().set_transform(Matrix::eye(4)).set_material(m1));

        let m2 = Material::new(Solid::new(Color::new(1., 1., 1.)), 1., 0.9, 0.9, 200.);
        let s2 = Rc::new(
            Sphere::new()
                .set_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .set_material(m2),
        );

        let w = World::new(vec![s1, s2], light);

        let inner_color = Color::new(1., 1., 1.);

        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let c = w.color_at(&r);
        assert_eq!(c, inner_color);
    }

    #[test]
    fn is_shadowed() {
        // there is no shadow when nothing is collinear with the point and light
        let w = World::default();
        let p = Point::new(0., 10., 0.);
        assert!(!w.is_shadowed(&p));

        // the shadow when an object is between the point and the light
        let p = Point::new(10., -10., 10.);
        assert!(w.is_shadowed(&p));

        // there is no shadow when an object is behind the light
        let p = Point::new(-20., 20., -20.);
        assert!(!w.is_shadowed(&p));

        // there is no shadow when an object is behind the point
        let p = Point::new(-2., 2., -2.);
        assert!(!w.is_shadowed(&p));
    }
}
