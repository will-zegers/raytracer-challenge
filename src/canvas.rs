use crate::color::Color;

pub struct Canvas {
    pub width: usize,
    pub height: usize,
    pixels: Vec<Color>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            pixels: vec![Color::new(0., 0., 0.); width * height],
        }
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, color: Color) {
        let idx = y * self.width + x;
        self.pixels[idx] = color;
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> &Color {
        let idx = y * self.width + x;
        &self.pixels[idx]
    }

    pub fn to_ppm(&self) -> String {
        let mut out: String = format!("P3\n{} {}\n255\n", self.width, self.height);
        for p in &self.pixels {
            let c = Color {
                r: f64::ceil(255. * Canvas::clamp(p.r)),
                g: f64::ceil(255. * Canvas::clamp(p.g)),
                b: f64::ceil(255. * Canvas::clamp(p.b)),
            };
            out += &format!("{} {} {}\n", c.r as u32, c.g as u32, c.b as u32);
        }
        out
    }

    fn clamp(n: f64) -> f64 {
        f64::min(1., f64::max(0., n))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new() {
        let c = Canvas::new(10, 20);
        assert_eq!(c.width, 10);
        assert_eq!(c.height, 20);
        for p in c.pixels {
            assert_eq!(p, Color::new(0., 0., 0.));
        }
    }

    #[test]
    fn write_pixel() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1., 0., 0.);
        c.write_pixel(2, 3, red);
        assert_eq!(*c.pixel_at(2, 3), Color::new(1., 0., 0.));
    }

    #[test]
    fn to_ppm() {
        let mut c = Canvas::new(5, 3);

        let c1 = Color::new(1.5, 0., 0.);
        let c2 = Color::new(0., 0.5, 0.);
        let c3 = Color::new(-0.5, 0., 1.);

        c.write_pixel(0, 0, c1);
        c.write_pixel(2, 1, c2);
        c.write_pixel(4, 2, c3);

        assert_eq!(
            c.to_ppm(),
            String::from(
                r#"P3
5 3
255
255 0 0
0 0 0
0 0 0
0 0 0
0 0 0
0 0 0
0 0 0
0 128 0
0 0 0
0 0 0
0 0 0
0 0 0
0 0 0
0 0 0
0 0 255
"#
            )
        );
    }
}
