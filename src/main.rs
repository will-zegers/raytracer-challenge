#![allow(unused_imports)]
use std::f64::consts::PI;

mod canvas;
use canvas::Canvas;

mod color;
use color::Color;

mod matrix;
use matrix::{Axis, Matrix};

mod point;
use point::Point;

mod vector;
use vector::Vector;

fn write_square(canvas: &mut Canvas, x: f64, y: f64) {
    let x = (x + (canvas.width as f64 / 2.)) as usize;
    let y = (y + (canvas.width as f64 / 2.)) as usize;
    let size = 3;

    for i in x-size .. x+size {
        for j in y-size .. y+size {
            canvas.write_pixel(i, j, Color::new(1., 1., 1.));
        }
    }
}

fn main() {
    let mut c = Canvas::new(600, 600);
    let mut angle = 0.;

    for _ in 0..12 {
        let t1 = Matrix::translation(0., -250., 0.);
        let t2 = Matrix::rotation(Axis::Z, angle);
        let p0 = Point::new(0., 0., 0.);

        let p1 = t2 * t1 * p0;
        write_square(&mut c, p1.x, p1.y);

        angle += PI / 6.;
    }

    println!("{}", c.to_ppm());
}
