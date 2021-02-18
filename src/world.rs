use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::color::Color;
use crate::hit_record::HitRecord;
use crate::intersection::IntersectionList;
use crate::light;
use crate::light::PointLight;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::ray::Ray;
use crate::sphere::Sphere;

#[allow(dead_code)]
pub struct World {
    objects: Vec<Sphere>,
    pub light: PointLight,
}

impl World {
    pub fn new(objects: Vec<Sphere>, light: PointLight) -> Self {
        Self {
            objects,
            light,
        }
    }

    #[allow(dead_code)]
    pub fn default() -> Self {
        let light = PointLight::default();
        let m1 = Material::new(Color::new(0.8, 1.0, 0.6), 0.1, 0.7, 0.2, 200.);
        let mut s1 = Sphere::default();
        s1.material = m1;
        let mut s2 = Sphere::default();
        s2.set_transform(Matrix::scaling(0.5, 0.5, 0.5));

        Self {
            objects: vec![s1, s2],
            light,
        }
    }

    #[allow(dead_code)]
    fn intersects(&self, r: &Ray) -> IntersectionList {
        let mut xs = IntersectionList::empty();
        for obj in &self.objects {
            match r.intersects(&obj) {
                Some(hits) => xs.extend(hits),
                None => (),
            }
        }
        xs
    }

    #[allow(dead_code)]
    fn shade(&self, rec: &HitRecord) -> Color {
        light::lighting(
            &rec.object.material,
            &self.light,
            &rec.point,
            &rec.eyev,
            &rec.normalv,
        )
    }

    #[allow(dead_code)]
    pub fn color_at(&self, r: &Ray) -> Color {
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
        for y in 0 .. c.vsize {
            for x in 0 .. c.hsize {
                let ray = c.pixel_to_ray(x, y);
                let color = self.color_at(&ray);
                image.write_pixel(x, y, color);
            }
        }
        image
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::intersection::Intersection;
    use crate::point::Point;
    use crate::vector::Vector;

    #[test]
    fn default() {
        let m1 = Material::new(Color::new(0.8, 1.0, 0.6), 0.1, 0.7, 0.2, 200.);
        let mut s1 = Sphere::default();
        s1.material = m1;
        let mut s2 = Sphere::default();
        s2.set_transform(Matrix::scaling(0.5, 0.5, 0.5));

        let w = World::default();
        assert_eq!(w.light, PointLight::default());
        assert!(w.objects.iter().any(|s| *s == s1));
        assert!(w.objects.iter().any(|s| *s == s2));
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
        assert_eq!(xs[0], Intersection::new(4.0, &w.objects[0]));

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
        let mut w = World::default();
        let inner_color;
        {
            let mut outer = w.objects.get_mut(0).unwrap();
            outer.material.ambient = 1.;
            let mut inner = w.objects.get_mut(1).unwrap();
            inner.material.ambient = 1.;
            inner_color = inner.material.color;
        }

        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let c = w.color_at(&r);
        assert_eq!(c, inner_color);
    }
}
