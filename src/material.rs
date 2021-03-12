use std::rc::Rc;

use crate::color::{Color, BLACK, WHITE};
use crate::patterns::{Pattern, Solid};

#[cfg_attr(test, derive(Debug))]
pub struct Material {
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub pattern: Box<dyn Pattern>,
    pub reflective: f64,
    pub transparency: f64,
    pub refractive_index: f64,
}

impl Material {
    pub fn new(
        pattern: impl Pattern + 'static,
        ambient: f64,
        diffuse: f64,
        specular: f64,
        shininess: f64,
        reflective: f64,
        transparency: f64,
        refractive_index: f64,
    ) -> Self {
        Self {
            pattern: Box::new(pattern),
            ambient,
            diffuse,
            specular,
            shininess,
            reflective,
            transparency,
            refractive_index,
        }
    }

    pub fn default() -> Self {
        Self {
            pattern: Box::new(Solid::new(Color::new(1., 1., 1.))),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.,
            reflective: 0.,
            transparency: 0.,
            refractive_index: 1.,
        }
    }

    pub fn with_pattern(mut self, pattern: impl Pattern + 'static) -> Self {
        self.pattern = Box::new(pattern);
        self
    }

    pub fn with_ambient(mut self, ambient: f64) -> Self {
        self.ambient = ambient;
        self
    }

    pub fn with_diffuse(mut self, diffuse: f64) -> Self {
        self.diffuse = diffuse;
        self
    }

    pub fn with_specular(mut self, specular: f64) -> Self {
        self.specular = specular;
        self
    }

    pub fn with_shininess(mut self, shininess: f64) -> Self {
        self.shininess = shininess;
        self
    }

    pub fn with_reflective(mut self, reflective: f64) -> Self {
        self.reflective = reflective;
        self
    }

    pub fn with_transparency(mut self, transparency: f64) -> Self {
        self.transparency = transparency;
        self
    }

    pub fn with_refractive_index(mut self, refractive_index: f64) -> Self {
        self.refractive_index = refractive_index;
        self
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use std::cmp::PartialEq;

    impl PartialEq for Material {
        fn eq(&self, other: &Self) -> bool {
            self.ambient == other.ambient
                && self.diffuse == other.diffuse
                && self.specular == other.specular
                && self.shininess == other.shininess
        }
    }

    #[test]
    fn new() {
        let m = Material::default();
        assert_eq!(m.ambient, 0.1);
        assert_eq!(m.diffuse, 0.9);
        assert_eq!(m.specular, 0.9);
        assert_eq!(m.shininess, 200.);
        assert_eq!(m.transparency, 0.);
        assert_eq!(m.refractive_index, 1.);
    }
}
