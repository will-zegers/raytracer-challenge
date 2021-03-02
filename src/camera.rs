use crate::matrix::Matrix;
use crate::point::Point;
use crate::ray::Ray;
use crate::vector::Vector;

type ViewTransform = Matrix;

pub fn view_transform(from: Point, to: Point, up: Vector) -> ViewTransform {
    let forward = (to - from).normalize();
    let upn = up.normalize();
    let left = Vector::cross(&forward, &upn);
    let true_up = Vector::cross(&left, &forward);

    let orientation = Matrix::new(
        4,
        4,
        // |     left.x |     left.y |     left.z | 0. |
        // |  true_up.x |  true_up.y |  true_up.z | 0. |
        // | -forward.x | -forward.y | -forward.z | 0. |
        // |        0.0 |        0.0 |        0.0 | 1. |
        vec![
            left.x, left.y, left.z, 0., true_up.x, true_up.y, true_up.z, 0., -forward.x,
            -forward.y, -forward.z, 0., 0., 0., 0., 1.,
        ],
    );

    orientation * Matrix::translation(-from.x, -from.y, -from.z)
}

pub struct Camera {
    pub hsize: usize,
    pub vsize: usize,
    pub pixel_size: f64,

    transform: Matrix,
    transform_inv: Matrix,
    half_width: f64,
    half_height: f64,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, fov: f64) -> Self {
        let aspect_ratio = (hsize as f64) / (vsize as f64);
        let half_view = f64::tan(fov / 2.);

        let (half_width, half_height): (f64, f64);
        if aspect_ratio >= 1. {
            half_width = half_view;
            half_height = half_view / aspect_ratio;
        } else {
            half_width = half_view * aspect_ratio;
            half_height = half_view;
        }
        let pixel_size = (half_width * 2.) / (hsize as f64);

        Self {
            hsize,
            vsize,
            pixel_size,
            transform: Matrix::eye(4),
            transform_inv: Matrix::eye(4),
            half_width,
            half_height,
        }
    }

    pub fn set_transform(mut self, t: Matrix) -> Self {
        self.transform_inv = t.inverse();
        self.transform = t;

        self
    }

    pub fn pixel_to_ray(&self, px: usize, py: usize) -> Ray {
        // the offset from the edge of the canvas to the pixel's center
        let xoffset = ((px as f64) + 0.5) * self.pixel_size;
        let yoffset = ((py as f64) + 0.5) * self.pixel_size;

        // the untransformed coordinates of the pixel in world space
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        // using the camera matrix, tranform the canvas point and the origin,
        // and then compute the ray's direction vector
        let pixel = &self.transform_inv * Point::new(world_x, world_y, -1.);
        let origin = &self.transform_inv * Point::new(0., 0., 0.);
        let direction = (pixel - origin).normalize();

        Ray::new(origin, direction)
    }
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use super::*;

    use crate::matrix::Axis;

    const TOL: f64 = 1e-9;

    fn is_close(a: f64, b: f64) -> bool {
        f64::abs(a - b) < TOL
    }

    #[test]
    fn view_transform() {
        // the trasformation matrix for the default orientation
        let from = Point::new(0., 0., 0.);
        let to = Point::new(0., 0., -1.);
        let up = Vector::new(0., 1., 0.);
        let t = super::view_transform(from, to, up);
        assert_eq!(t, Matrix::eye(4));

        // a ViewTransform looking in the positive z direction
        let from = Point::new(0., 0., 0.);
        let to = Point::new(0., 0., 1.);
        let up = Vector::new(0., 1., 0.);
        let t = super::view_transform(from, to, up);
        assert_eq!(t, Matrix::scaling(-1., 1., -1.));

        // the view transformation moves the world
        let from = Point::new(0., 0., 8.);
        let to = Point::new(0., 0., 0.);
        let up = Vector::new(0., 1., 0.);
        let t = super::view_transform(from, to, up);
        assert_eq!(t, Matrix::translation(0., 0., -8.));

        // an arbitrary view transformation
        let from = Point::new(1., 3., 2.);
        let to = Point::new(4., -2., 8.);
        let up = Vector::new(1., 1., 0.);
        let t = super::view_transform(from, to, up);
        assert_eq!(
            t,
            Matrix::new(
                4,
                4,
                vec![
                    -0.507092552,
                    0.507092552,
                    0.676123403,
                    -2.366431913,
                    0.767715933,
                    0.606091526,
                    0.121218305,
                    -2.828427124,
                    -0.358568582,
                    0.597614304,
                    -0.717137165,
                    0.,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                ]
            )
        )
    }

    #[test]
    fn new() {
        let hsize = 160;
        let vsize = 120;
        let fov = PI / 2.;
        let c = Camera::new(hsize, vsize, fov);

        assert_eq!(c.hsize, hsize);
        assert_eq!(c.vsize, vsize);
        assert_eq!(c.transform, Matrix::eye(4));
    }

    #[test]
    fn pixel_size() {
        let c = Camera::new(200, 125, PI / 2.);
        assert!(is_close(c.pixel_size, 0.01), "{} != 0.01", c.pixel_size);

        let c = Camera::new(125, 200, PI / 2.);
        assert!(is_close(c.pixel_size, 0.01), "{} != 0.01", c.pixel_size);
    }

    #[test]
    pub fn pixel_to_ray() {
        // constructing a ray through the center of a canvas
        let c = Camera::new(201, 101, PI / 2.);
        let r = c.pixel_to_ray(100, 50);
        assert_eq!(r.origin, Point::new(0., 0., 0.));
        assert_eq!(r.direction, Vector::new(0., 0., -1.));

        // constructing a ray through the corner of a canvas
        let c = Camera::new(201, 101, PI / 2.);
        let r = c.pixel_to_ray(0, 0);
        assert_eq!(r.origin, Point::new(0., 0., 0.));
        assert_eq!(
            r.direction,
            Vector::new(0.665186426, 0.332593213, -0.668512358)
        );

        // constructing a ray when the camera is transformed
        let c = Camera::new(201, 101, PI / 2.)
            .set_transform(Matrix::rotation(Axis::Y, PI / 4.) * Matrix::translation(0., -2., 5.));
        let r = c.pixel_to_ray(100, 50);
        assert_eq!(r.origin, Point::new(0., 2., -5.));
        assert_eq!(
            r.direction,
            Vector::new(f64::sqrt(2.) / 2., 0., -f64::sqrt(2.) / 2.)
        );
    }
}
