use std::rc::Rc;

use crate::color::Color;
use crate::material::Material;
use crate::patterns;
use crate::patterns::Pattern;
use crate::point::Point;
use crate::shape::Shape;
use crate::vector::Vector;

#[cfg_attr(test, derive(Debug, PartialEq))]
pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            position,
            intensity,
        }
    }

    pub fn default() -> Self {
        Self {
            position: Point::new(-10., 10., -10.),
            intensity: Color::new(1., 1., 1.),
        }
    }
}

pub fn lighting(
    obj: Rc<dyn Shape>,
    light: &PointLight,
    point: &Point,
    eyev: &Vector,
    normalv: &Vector,
    in_shadow: bool,
) -> Color {
    // combine the surface color with the light's color/intensity
    let material = obj.material();
    let effective_color = material.pattern.color_at_object(obj.clone(), &point) * light.intensity;

    // find the direction to the light source
    let lightv = (light.position - *point).normalize();

    // compute the ambient contribution
    let ambient = effective_color * material.ambient;

    let mut diffuse = Color::new(0., 0., 0.);
    let mut specular = Color::new(0., 0., 0.);

    if !in_shadow {
        // light_dot_normal represent the cosine of the angle between the light
        // vector and the normal. A negative number means the light is on the other
        // side of the surface
        let light_dot_normal = Vector::dot(&lightv, normalv);
        if light_dot_normal >= 0. {
            diffuse = effective_color * material.diffuse * light_dot_normal;

            // reflect_dot_eye represents the cosine of the angle between the
            // reflection vector and the eye vector. A negative number means the
            // light reflects away from the eye
            let reflectv = -lightv.reflect(normalv);
            let reflect_dot_eye = Vector::dot(&reflectv, &eyev.normalize());
            if reflect_dot_eye > 0. {
                let factor = f64::powf(reflect_dot_eye, material.shininess);
                specular = light.intensity * material.specular * factor;
            }
        }
    }
    ambient + diffuse + specular
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::color::{BLACK, WHITE};
    use crate::patterns::Stripe;
    use crate::shape::Sphere;

    #[test]
    fn new() {
        let position = Point::new(0., 0., 0.);
        let intensity = WHITE;
        let light = PointLight::new(position, intensity);
        assert_eq!(light.position, position);
        assert_eq!(light.intensity, intensity);
    }

    #[test]
    fn default() {
        let light = PointLight::default();
        assert_eq!(light.position, Point::new(10., 10., 10.));
        assert_eq!(light.intensity, WHITE);
    }

    #[test]
    fn test_lighting() {
        let m = Material::default();
        let sphere = Rc::new(Sphere::new().with_material(Material::default()));
        let position = Point::new(0., 0., 0.);

        // lighting with the eye between the light and the surface
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), WHITE);
        let result = lighting(sphere.clone(), &light, &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));

        // light with the eye between light and surface, eye offset 45⁰
        let eyev = Vector::new(0., f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.);
        let result = lighting(sphere.clone(), &light, &position, &eyev, &normalv, false);
        assert_eq!(result, WHITE);

        // light with the eye between light and surface, light offset 45⁰
        let eyev = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), WHITE);
        let result = lighting(sphere.clone(), &light, &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(0.736396103, 0.736396103, 0.736396103));

        // lighting with eye in the path of the reflection vector
        let eyev = Vector::new(0., -f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.);
        let result = lighting(sphere.clone(), &light, &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(1.636396103, 1.636396103, 1.636396103));

        // lighting with the light behind the surface
        let eyev = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., 10.), WHITE);
        let result = lighting(sphere, &light, &position, &eyev, &normalv, false);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_shadow() {
        let m = Material::default();
        let sphere = Rc::new(Sphere::new().with_material(Material::default()));
        let position = Point::new(0., 0., 0.);

        // lighting with the surface in shadow
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), WHITE);
        let in_shadow = true;

        let result = super::lighting(sphere, &light, &position, &eyev, &normalv, in_shadow);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn light_with_pattern() {
        let m = Material::new(Stripe::new(WHITE, BLACK), 1., 0., 0., 200., 0., 0., 1.);
        let sphere = Rc::new(Sphere::new().with_material(m));
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), WHITE);
        let c1 = lighting(
            sphere.clone(),
            &light,
            &Point::new(0.9, 0., 0.),
            &eyev,
            &normalv,
            false,
        );
        let c2 = lighting(
            sphere.clone(),
            &light,
            &Point::new(1.1, 0., 0.),
            &eyev,
            &normalv,
            false,
        );
        assert_eq!(c1, WHITE);
        assert_eq!(c2, BLACK);
    }
}
