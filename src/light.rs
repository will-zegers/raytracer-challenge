use crate::color::Color;
use crate::material::Material;
use crate::point::Point;
use crate::vector::Vector;

#[derive(Debug, PartialEq)]
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

    #[allow(dead_code)]
    pub fn default() -> Self {
        Self {
            position: Point::new(-10., 10., -10.),
            intensity: Color::new(1., 1., 1.),
        }
    }
}

pub fn lighting(
    m: &Material,
    light: &PointLight,
    point: &Point,
    eyev: &Vector,
    normalv: &Vector,
) -> Color {
    // combine the surface color with the light's color/intensity
    let effective_color = m.color * light.intensity;

    // find the direction to the light source
    let lightv = (light.position - *point).normalize();

    // compute the ambient contribution
    let ambient = effective_color * m.ambient;

    let mut diffuse = Color::new(0., 0., 0.);
    let mut specular = Color::new(0., 0., 0.);

    // light_dot_normal represent the cosine of the angle between the light
    // vector and the normal. A negative number means the light is on the other
    // side of the surface
    let light_dot_normal = Vector::dot(&lightv, normalv);
    if light_dot_normal >= 0. {
        diffuse = effective_color * m.diffuse * light_dot_normal;

        // reflect_dot_eye represents the cosine of the angle between the
        // reflection vector and the eye vector. A negative number means the
        // light reflects away from the eye
        let reflectv = -lightv.reflect(normalv);
        let reflect_dot_eye = Vector::dot(&reflectv, &eyev.normalize());
        if reflect_dot_eye > 0. {
            let factor = f64::powf(reflect_dot_eye, m.shininess);
            specular = light.intensity * m.specular * factor;
        }
    }
    ambient + diffuse + specular
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new() {
        let position = Point::new(0., 0., 0.);
        let intensity = Color::new(1., 1., 1.);
        let light = PointLight::new(position, intensity);
        assert_eq!(light.position, position);
        assert_eq!(light.intensity, intensity);
    }

    #[test]
    fn default() {
        let light = PointLight::default();
        assert_eq!(light.position, Point::new(10., 10., 10.));
        assert_eq!(light.intensity, Color::new(1., 1., 1.));
    }

    #[test]
    fn test_lighting() {
        let m = Material::default();
        let position = Point::new(0., 0., 0.);

        // lighting with the eye between the light and the surface
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::new(1., 1., 1.));
        let result = lighting(&m, &light, &position, &eyev, &normalv);
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));

        // light with the eye between light and surface, eye offset 45⁰
        let eyev = Vector::new(0., f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.);
        let result = lighting(&m, &light, &position, &eyev, &normalv);
        assert_eq!(result, Color::new(1.0, 1.0, 1.0));

        // light with the eye between light and surface, light offset 45⁰
        let eyev = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::new(1., 1., 1.));
        let result = lighting(&m, &light, &position, &eyev, &normalv);
        assert_eq!(result, Color::new(0.736396103, 0.736396103, 0.736396103));

        // lighting with eye in the path of the reflection vector
        let eyev = Vector::new(0., -f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.);
        let result = lighting(&m, &light, &position, &eyev, &normalv);
        assert_eq!(result, Color::new(1.636396103, 1.636396103, 1.636396103));

        // lighting with the light behind the surface
        let eyev = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., 10.), Color::new(1., 1., 1.));
        let result = lighting(&m, &light, &position, &eyev, &normalv);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }
}
