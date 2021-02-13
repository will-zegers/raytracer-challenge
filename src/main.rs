#![allow(unused_imports)]
#![allow(unused_variables)]
use std::f64::consts::PI;

mod canvas;
use canvas::Canvas;

mod color;
use color::Color;

mod intersection;

mod light;
use light::PointLight;

mod material;

mod matrix;
use matrix::{Axis, Matrix};

mod point;
use point::Point;

mod ray;
use crate::ray::Ray;

mod sphere;
use crate::sphere::Sphere;

mod vector;
use vector::Vector;

fn main() {
    let (width, height) = (600, 600);
    let aspect_ratio = (width as f64) / (height as f64);
    let mut cv = Canvas::new(width, height);

    let mut s = Sphere::new();
    s.material.color = Color::new(1., 0.2, 1.);

    let light = PointLight::new(Point::new(-10., 10., -10.), Color::new(1., 1., 1.));

    let origin = Point::new(0., 0., -5.);
    let z_wall = if origin.z < 0. {
        -origin.z + 1.
    } else {
        -origin.z - 1.
    };
    let wall_width = 4.;
    for i in 0..height {
        let y = (wall_width / (2. * aspect_ratio)) - ((i as f64) * wall_width) / (aspect_ratio * (height as f64));

        for j in 0..width {
            let x = ((j as f64) * wall_width) / (width as f64) - (wall_width / 2.);

            let r = Ray::new(origin, Vector::new(x, y, z_wall));
            match r.intersects(&s) {
                Some(hit) => {
                    let hit = &hit[0];
                    let point = r.at(hit.t);
                    let normal = hit.object.normal_at(point);
                    let eye = -r.direction;

                    let color =
                        light::lighting(&hit.object.material, &light, &point, &eye, &normal);
                    cv.write_pixel(j, i, color);
                }
                None => (),
            }
        }
    }
    println!("{}", cv.to_ppm());
}
