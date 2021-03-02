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
}

impl Material {
    pub fn new(
        pattern: impl Pattern + 'static,
        ambient: f64,
        diffuse: f64,
        specular: f64,
        shininess: f64,
    ) -> Self {
        Self {
            pattern: Box::new(pattern),
            ambient,
            diffuse,
            specular,
            shininess,
        }
    }

    pub fn default() -> Self {
        Self {
            pattern: Box::new(Solid::new(Color::new(1., 1., 1.))),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.,
        }
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
    }
}
