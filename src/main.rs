#![allow(unused_imports)]
#![allow(unused_variables)]
use std::f64::consts::PI;
use std::rc::Rc;

mod camera;
use camera::Camera;

mod canvas;
use canvas::Canvas;

mod color;
use color::{Color, BLACK, WHITE};

mod hit_record;

mod intersection;

mod light;
use light::PointLight;

mod material;
use material::Material;

mod matrix;
use matrix::{Axis, Matrix};

mod patterns;
use patterns::{Checker, Gradient, Pattern, Ring, Solid, Stripe};

mod point;
use point::Point;

mod ray;
use ray::Ray;

mod shape;
use shape::{Plane, Shape, Sphere};

mod vector;
use vector::Vector;

mod world;
use world::World;

fn main() {
    let mut objects = Vec::<Rc<dyn Shape>>::new();

    let floor_mat = Material::new(Checker::new(WHITE, BLACK), 0.1, 0.9, 0., 200., 0.5, 0., 1.);
    let floor = Plane::new(Matrix::eye(4), floor_mat);
    objects.push(Rc::new(floor));

    let wall_mat = Material::new(
        Solid::new(Color::new(0.5, 0.7, 0.9)),
        0.0,
        0.5,
        0.,
        0.,
        0.,
        0.,
        1.,
    );
    let wall_tf = Matrix::translation(0., 0., 15.) * Matrix::rotation(Axis::X, PI / 2.);
    let wall = Plane::new(wall_tf, wall_mat);
    objects.push(Rc::new(wall));

    let middle_tf = Matrix::translation(-0.5, 1., 0.5);
    let middle_mat = Material::new(
        Solid::new(Color::new(0.6, 0.6, 0.6)),
        0.1,
        0.7,
        0.0,
        200.,
        0.2,
        0.,
        1.,
    );
    let middle = Sphere::new()
        .set_transform(middle_tf)
        .set_material(middle_mat);
    objects.push(Rc::new(middle));

    let right_tf = Matrix::translation(1.75, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    let right_pattern = Checker::new(Color::new(0.9, 0.1, 0.1), Color::new(0.1, 0.1, 0.9))
        .set_transform(Matrix::rotation(Axis::Z, PI / 4.) * Matrix::scaling(0.3, 0.3, 0.3));
    let right_mat = Material::new(right_pattern, 0.1, 0.7, 0.3, 100., 0., 0., 1.);
    let right = Sphere::new()
        .set_transform(right_tf)
        .set_material(right_mat);
    objects.push(Rc::new(right));

    let left_tf = Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33);
    let left_pattern = Ring::new(WHITE, BLACK)
        .set_transform(Matrix::rotation(Axis::X, PI / 2.) * Matrix::scaling(0.1, 0.1, 0.1));
    let left_mat = Material::new(left_pattern, 0.1, 0.7, 0.3, 300., 0., 0., 1.);
    let left = Sphere::new().set_transform(left_tf).set_material(left_mat);
    objects.push(Rc::new(left));

    let camera_tf = camera::view_transform(
        Point::new(0., 1.5, -5.),
        Point::new(0., 1., 0.),
        Vector::new(0., 1., 0.),
    );
    let camera = Camera::new(848, 480, PI / 3.).set_transform(camera_tf);

    let world = World::new(objects, PointLight::default());
    let canvas = world.render(camera, 5);
    println!("{}", canvas.to_ppm());
}
