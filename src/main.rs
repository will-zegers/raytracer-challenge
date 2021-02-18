#![allow(unused_imports)]
#![allow(unused_variables)]
use std::f64::consts::PI;

mod camera;
use camera::Camera;

mod canvas;
use canvas::Canvas;

mod color;
use color::Color;

mod hit_record;

mod intersection;

mod light;
use light::PointLight;

mod material;
use material::Material;

mod matrix;
use matrix::{Axis, Matrix};

mod point;
use point::Point;

mod ray;
use ray::Ray;

mod sphere;
use sphere::Sphere;

mod vector;
use vector::Vector;

mod world;
use world::World;

fn main() {
    let mut objects = Vec::<Sphere>::new();

    let floor_mat = Material::new(Color::new(1., 0.9, 0.9), 0.1, 0.9, 0., 200.);
    let floor = Sphere::new(Matrix::scaling(10., 0.01, 10.), floor_mat);
    objects.push(floor);

    let left_wall_tf = Matrix::translation(0., 0., 5.)
        * Matrix::rotation(Axis::Y, -PI / 4.)
        * Matrix::rotation(Axis::X, PI / 2.)
        * Matrix::scaling(10., 0.01, 10.);
    let left_wall = Sphere::new(left_wall_tf, floor_mat);
    objects.push(left_wall);

    let right_wall_tf = Matrix::translation(0., 0., 5.)
        * Matrix::rotation(Axis::Y, PI / 4.)
        * Matrix::rotation(Axis::X, PI / 2.)
        * Matrix::scaling(10., 0.01, 10.);
    let right_wall = Sphere::new(right_wall_tf, floor_mat);
    objects.push(right_wall);

    let middle_tf = Matrix::translation(-0.5, 1., 0.5);
    let middle_mat = Material::new(Color::new(0.1, 1., 0.5), 0.1, 0.7, 0.3, 200.);
    let middle = Sphere::new(middle_tf, middle_mat);
    objects.push(middle);

    let right_tf = Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    let right_mat = Material::new(Color::new(0.5, 1., 0.1), 0.1, 0.7, 0.3, 200.);
    let right = Sphere::new(right_tf, right_mat);
    objects.push(right);

    let left_tf = Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33);
    let left_mat = Material::new(Color::new(1., 0.8, 0.1), 0.1, 0.7, 0.3, 200.);
    let left = Sphere::new(left_tf, left_mat);
    objects.push(left);

    let mut camera = Camera::new(1920, 1080, PI / 3.);
    let camera_tf = camera::view_transform(Point::new(0., 1.5, -5.), Point::new(0., 1., 0.), Vector::new(0., 1., 0.));
    camera.set_transform(camera_tf);

    let world = World::new(objects, PointLight::default());
    let canvas = world.render(camera);
    println!("{}", canvas.to_ppm());
}

