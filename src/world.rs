use std::rc::Rc;

use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::color;
use crate::color::Color;
use crate::hit_record::HitRecord;
use crate::light;
use crate::light::PointLight;
use crate::material::Material;
use crate::matrix::Matrix;
use crate::patterns::Solid;
use crate::point::Point;
use crate::ray::Ray;
use crate::shape::{Shape, Sphere};
use crate::vector::Vector;
use crate::{
    intersection,
    intersection::{Intersection, IntersectionList},
};

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
        let m1 = Material::new(
            Solid::new(Color::new(0.8, 1.0, 0.6)),
            0.1,
            0.7,
            0.2,
            200.,
            0.,
            0.,
            1.,
        );
        let s1 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::eye(4))
                .with_material(m1),
        );
        let s2 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .with_material(Material::default()),
        );

        Self {
            objects: vec![s1, s2],
            light,
        }
    }

    pub fn add(&mut self, shape: Rc<dyn Shape>) {
        self.objects.push(shape);
    }

    pub fn intersects(&self, r: Ray) -> IntersectionList {
        let mut xs = IntersectionList::new();
        for obj in &self.objects {
            let hits = r.intersects(obj.clone());
            if !hits.is_empty() {
                xs.extend(hits);
            }
        }
        intersection::sort(&mut xs);
        xs
    }

    fn shade(&self, rec: &HitRecord, remaining: i32) -> Color {
        let surface = light::lighting(
            rec.object.clone(),
            &self.light,
            &rec.over_point,
            &rec.eyev,
            &rec.normalv,
            self.is_shadowed(&rec.over_point),
        );
        let reflected = self.reflected_color(rec, remaining);
        let refracted = self.refracted_color(rec, remaining);

        if rec.object.material().reflective > 0. && rec.object.material().transparency > 0. {
            let reflectance = rec.schlick();
            return surface + reflected * reflectance + refracted * (1. - reflectance);
        }

        surface + reflected + refracted
    }

    fn color_at(&self, r: Ray, remaining: i32) -> Color {
        let mut xs = self.intersects(r);
        intersection::sort(&mut xs);
        return match intersection::hit(&xs) {
            Some(x) => {
                let rec = HitRecord::new(x, r, &xs);
                self.shade(&rec, remaining)
            }
            None => Color::new(0., 0., 0.),
        };
    }

    pub fn render(&self, c: Camera, depth: i32) -> Canvas {
        let mut image = Canvas::new(c.hsize, c.vsize);
        for y in 0..c.vsize {
            for x in 0..c.hsize {
                let ray = c.pixel_to_ray(x, y);
                let color = self.color_at(ray, depth);
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
        let mut xs = self.intersects(r);

        return match intersection::hit(&mut xs) {
            Some(hit) => hit.t < distance,
            None => false,
        };
    }

    fn reflected_color(&self, rec: &HitRecord, remaining: i32) -> Color {
        if remaining == 0 {
            return color::BLACK;
        }
        let r = Ray::new(rec.over_point, rec.reflectv);
        self.color_at(r, remaining - 1) * rec.object.material().reflective
    }

    fn refracted_color(&self, rec: &HitRecord, remaining: i32) -> Color {
        if remaining == 0 || rec.object.material().transparency == 0. {
            return color::BLACK;
        }

        // create the refracted ray
        return match Vector::refract(&rec.normalv, &rec.eyev, rec.n1 / rec.n2) {
            Some(direction) => {
                let refract_ray = Ray::new(rec.under_point, direction);
                self.color_at(refract_ray, remaining - 1) * rec.object.material().transparency
            }
            None => color::BLACK,
        };
    }
}

#[cfg(test)]
mod test {
    use std::rc::Rc;

    use super::*;

    use crate::intersection::Intersection;
    use crate::patterns::TestPattern;
    use crate::point::Point;
    use crate::shape::Plane;
    use crate::vector::Vector;

    const DEFAULT_DEPTH: i32 = 5;

    #[test]
    fn default() {
        let m1 = Material::new(
            Solid::new(Color::new(0.8, 1.0, 0.6)),
            0.1,
            0.7,
            0.2,
            200.,
            0.,
            0.,
            1.,
        );
        let s1 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::eye(4))
                .with_material(m1),
        );
        let s2 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .with_material(Material::default()),
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
        let xs = w.intersects(r);
        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].t, 4.);
        assert_eq!(xs[1].t, 4.5);
        assert_eq!(xs[2].t, 5.5);
        assert_eq!(xs[3].t, 6.);
    }

    mod shade {
        use super::*;

        #[test]
        fn at_an_intersection() {
            let w = World::default();
            let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
            let xs = w.intersects(r);
            assert_ne!(xs.len(), 0);
            assert_eq!(xs[0], Intersection::new(4.0, w.objects[0].clone()));

            let rec = HitRecord::new(&xs[0], r, &mut IntersectionList::new());
            let c = w.shade(&rec, DEFAULT_DEPTH);
            assert_eq!(c, Color::new(0.380661193, 0.475826491, 0.285495894));
        }

        #[test]
        fn at_an_intersection_from_the_inside() {
            let mut w = World::default();
            w.light = PointLight::new(Point::new(0., 0.25, 0.), Color::new(1., 1., 1.));
            let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
            let xs = w.intersects(r);
            let rec = HitRecord::new(&xs[2], r, &mut IntersectionList::new());
            let c = w.shade(&rec, DEFAULT_DEPTH);
            assert_eq!(c, Color::new(0.904984472, 0.904984472, 0.904984472));
        }

        #[test]
        fn given_an_intersection_in_shadow() {
            let s1 = Rc::new(Sphere::new());
            let t2 = Matrix::translation(0., 0., 10.);
            let s2 = Rc::new(
                Sphere::new()
                    .with_transform(t2)
                    .with_material(Material::default()),
            );
            let light = PointLight::new(Point::new(0., 0., -10.), Color::new(1., 1., 1.));
            let w = World::new(vec![s1, s2], light);

            let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
            let mut xs = w.intersects(r);
            let hit = intersection::hit(&mut xs).unwrap();
            let rec = HitRecord::new(&hit, r, &mut IntersectionList::new());
            let c = w.shade(&rec, DEFAULT_DEPTH);
            assert_eq!(c, Color::new(0.1, 0.1, 0.1));
        }

        #[test]
        fn reflective_material() {
            let light = PointLight::default();

            let m1 = Material::new(
                Solid::new(Color::new(0.8, 1.0, 0.6)),
                0.1,
                0.7,
                0.2,
                200.,
                0.,
                0.,
                1.,
            );
            let s1 = Rc::new(
                Sphere::new()
                    .with_transform(Matrix::eye(4))
                    .with_material(m1),
            );

            let s2 = Rc::new(
                Sphere::new()
                    .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                    .with_material(Material::default()),
            );

            let mut mat_p1 = Material::default();
            mat_p1.reflective = 0.5;
            let p1 = Rc::new(
                Plane::default()
                    .with_transform(Matrix::translation(0., -1., 0.))
                    .with_material(mat_p1),
            );

            let w = World::new(vec![s1, s2, p1], light);

            let r = Ray::new(
                Point::new(0., 0., 3.),
                Vector::new(0., -f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.),
            );
            let mut xs = w.intersects(r);
            let i = intersection::hit(&mut xs).unwrap();
            let rec = HitRecord::new(&i, r, &mut IntersectionList::new());
            let color = w.shade(&rec, DEFAULT_DEPTH);
            assert_eq!(color, Color::new(0.584805085, 0.584805085, 0.584805085));
        }

        #[test]
        fn transparent_material() {
            let mut w = World::default();
            let mut m_floor = Material::default();
            m_floor.transparency = 0.5;
            m_floor.refractive_index = 1.5;
            let floor = Plane::default()
                .with_transform(Matrix::translation(0., -1., 0.))
                .with_material(m_floor);
            w.add(Rc::new(floor));

            let mut m_ball = Material::default();
            m_ball.pattern = Box::new(Solid::new(Color::new(1., 0., 0.)));
            m_ball.ambient = 0.5;
            let ball = Sphere::new()
                .with_transform(Matrix::translation(0., -3.5, -0.5))
                .with_material(m_ball);
            w.add(Rc::new(ball));

            let r = Ray::new(
                Point::new(0., 0., -3.),
                Vector::new(0., -f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.),
            );
            let xs = w.intersects(r);
            let rec = HitRecord::new(&xs[0], r, &xs);
            assert_eq!(
                w.shade(&rec, 5),
                Color::new(0.936425388, 0.686425388, 0.686425388)
            );
        }

        #[test]
        fn reflective_transparent_material() {
            let mut w = World::default();
            let r = Ray::new(
                Point::new(0., 0., -3.),
                Vector::new(0., -2_f64.sqrt() / 2., 2_f64.sqrt() / 2.),
            );

            let mut m_floor = Material::default();
            m_floor.reflective = 0.5;
            m_floor.transparency = 0.5;
            m_floor.refractive_index = 1.5;
            let floor = Plane::default()
                .with_material(m_floor)
                .with_transform(Matrix::translation(0., -1., 0.));
            w.add(Rc::new(floor));

            let mut m_ball = Material::default();
            m_ball.pattern = Box::new(Solid::new(Color::new(1., 0., 0.)));
            m_ball.ambient = 0.5;
            let ball = Sphere::new()
                .with_transform(Matrix::translation(0., -3.5, -0.5))
                .with_material(m_ball);
            w.add(Rc::new(ball));

            let xs = w.intersects(r);
            let rec = HitRecord::new(&xs[0], r, &xs);
            assert_eq!(
                w.shade(&rec, 5),
                Color::new(0.933915140, 0.696434226, 0.692430691)
            );
        }
    }

    #[test]
    fn color_at() {
        // the color when a ray misses
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        let c = w.color_at(r, DEFAULT_DEPTH);
        assert_eq!(c, Color::new(0., 0., 0.));

        // the color when a ray hits
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let c = w.color_at(r, DEFAULT_DEPTH);
        assert_eq!(c, Color::new(0.380661193, 0.475826491, 0.285495894));

        // the color with an intersection behind the ray
        let light = PointLight::default();

        let m1 = Material::new(
            Solid::new(Color::new(0.8, 1.0, 0.6)),
            1.,
            0.7,
            0.2,
            200.,
            0.,
            0.,
            1.,
        );
        let s1 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::eye(4))
                .with_material(m1),
        );

        let m2 = Material::new(
            Solid::new(Color::new(1., 1., 1.)),
            1.,
            0.9,
            0.9,
            200.,
            0.,
            0.,
            1.,
        );
        let s2 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .with_material(m2),
        );

        let w = World::new(vec![s1, s2], light);

        let inner_color = Color::new(1., 1., 1.);

        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let c = w.color_at(r, DEFAULT_DEPTH);
        assert_eq!(c, inner_color);

        // color_at with mutually reflected surfaces
        let mut mat_p1 = Material::default();
        mat_p1.reflective = 1.;
        let p1 = Rc::new(
            Plane::default()
                .with_transform(Matrix::translation(0., -1., 0.))
                .with_material(mat_p1),
        );

        let mut mat_p2 = Material::default();
        mat_p2.reflective = 1.;
        let p2 = Rc::new(
            Plane::default()
                .with_transform(Matrix::translation(0., 1., 0.))
                .with_material(mat_p2),
        );

        let light = PointLight::default();
        let w = World::new(vec![p1, p2], light);

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        w.color_at(r, DEFAULT_DEPTH); // should terminate, i.e. not recurse infinitely
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

    #[test]
    fn reflected_color() {
        let light = PointLight::default();
        let m1 = Material::new(
            Solid::new(Color::new(0.8, 1.0, 0.6)),
            0.1,
            0.7,
            0.2,
            200.,
            0.,
            0.,
            1.,
        );
        let s1 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::eye(4))
                .with_material(m1),
        );
        let m2 = Material::new(
            Solid::new(Color::new(1., 1., 1.)),
            1.,
            0.9,
            0.9,
            200.,
            0.,
            0.,
            1.,
        );
        let s2 = Rc::new(
            Sphere::new()
                .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                .with_material(Material::default()),
        );
        let mut w = World::new(vec![s1, s2], light);

        // strike a non-reflective surface
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(1., w.objects[1].clone());
        let rec = HitRecord::new(&i, r, &mut IntersectionList::new());
        let color = w.reflected_color(&rec, DEFAULT_DEPTH);

        assert_eq!(color, Color::new(0., 0., 0.));

        // strike a reflective surface
        let mut mat_p1 = Material::default();
        mat_p1.reflective = 0.5;
        let p1 = Plane::default()
            .with_transform(Matrix::translation(0., -1., 0.))
            .with_material(mat_p1);

        w.add(Rc::new(p1));

        let r = Ray::new(
            Point::new(0., 0., -3.),
            Vector::new(0., -f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.),
        );

        let mut xs = w.intersects(r);
        let i = intersection::hit(&mut xs).unwrap();
        let rec = HitRecord::new(&i, r, &mut IntersectionList::new());
        let color = w.reflected_color(&rec, DEFAULT_DEPTH);
        assert_eq!(color, Color::new(0.190330596, 0.237913245, 0.142747947));

        // the reflective color at maximum recusive depth
        let mut mat_p1 = Material::default();
        mat_p1.reflective = 0.5;
        let p1 = Plane::default()
            .with_transform(Matrix::translation(0., -1., 0.))
            .with_material(mat_p1);

        w.add(Rc::new(p1));

        let r = Ray::new(
            Point::new(0., 0., -3.),
            Vector::new(0., -f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.),
        );

        let mut xs = w.intersects(r);
        let i = intersection::hit(&mut xs).unwrap();
        let rec = HitRecord::new(&i, r, &mut IntersectionList::new());
        let color = w.reflected_color(&rec, 0);
        assert_eq!(color, Color::new(0., 0., 0.));
    }

    mod refracted_color {

        use super::*;

        #[test]
        fn with_an_opaque_surface() {
            let w = World::default();
            let s = w.objects.first().unwrap();
            let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
            let xs = w.intersects(r);

            let rec = HitRecord::new(&xs[0], r, &xs);
            assert_eq!(w.refracted_color(&rec, 5), color::BLACK);
        }

        #[test]
        fn at_maximum_recursive_depth() {
            let mut m1 = Material::default();
            m1.diffuse = 0.7;
            m1.specular = 0.2;
            m1.transparency = 1.;
            m1.refractive_index = 1.5;
            let s1 = Rc::new(
                Sphere::new()
                    .with_transform(Matrix::eye(4))
                    .with_material(m1),
            );

            let s2 = Rc::new(
                Sphere::new()
                    .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                    .with_material(Material::default()),
            );
            let w = World::new(vec![s1, s2], PointLight::default());

            let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
            let shape = w.objects.first().unwrap();
            let xs = w.intersects(r);
            let rec = HitRecord::new(&xs[0], r, &xs);
            assert_eq!(w.refracted_color(&rec, 0), color::BLACK);
        }

        #[test]
        fn under_total_internal_reflection() {
            let mut m1 = Material::default();
            m1.diffuse = 0.7;
            m1.specular = 0.2;
            m1.transparency = 1.;
            m1.refractive_index = 1.5;
            let s1 = Rc::new(
                Sphere::new()
                    .with_transform(Matrix::eye(4))
                    .with_material(m1),
            );

            let s2 = Rc::new(
                Sphere::new()
                    .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
                    .with_material(Material::default()),
            );
            let w = World::new(vec![s1, s2], PointLight::default());

            let r = Ray::new(
                Point::new(0., 0., f64::sqrt(2.) / 2.),
                Vector::new(0., 1., 0.),
            );
            let xs = w.intersects(r);
            let rec = HitRecord::new(&xs[1], r, &xs);
            assert_eq!(w.refracted_color(&rec, 5), color::BLACK);
        }

        #[test]
        fn refracted_ray() {
            let tf1 = Matrix::eye(4);
            let mut m1 = Material::default();
            m1.pattern = Box::new(TestPattern::new());
            m1.ambient = 1.;
            m1.specular = 0.2;
            let s1 = Rc::new(Sphere::new().with_transform(tf1).with_material(m1));

            let tf2 = Matrix::scaling(0.5, 0.5, 0.5);
            let mut m2 = Material::default();
            m2.transparency = 1.;
            m2.refractive_index = 1.5;
            let s2 = Rc::new(Sphere::new().with_transform(tf2).with_material(m2));
            let w = World::new(vec![s1.clone(), s2.clone()], PointLight::default());

            let r = Ray::new(Point::new(0., 0., 0.1), Vector::new(0., 1., 0.));
            let xs = w.intersects(r);
            let rec = HitRecord::new(&xs[2], r, &xs);
            assert_eq!(
                w.refracted_color(&rec, 5),
                Color::new(0.0, 0.998884703, 0.047215978)
            );
        }
    }
}
