#![allow(unused_imports)]
#![allow(unused_variables)]
use std::f64::consts::PI;
use std::rc::Rc;

mod camera;
use camera::Camera;

mod canvas;
use canvas::Canvas;

mod color;
use color::{Color, BLACK, DBROWN, LBROWN, WHITE};

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
use shape::{Cube, Plane, Shape, Sphere};

mod vector;
use vector::Vector;

mod world;
use world::World;

fn main() {
    let mut objects = Vec::<Rc<dyn Shape>>::new();

    let room_mat = Material::default()
        .with_pattern(Checker::new(BLACK, WHITE).with_transform(Matrix::scaling(0.2, 0.2, 0.2)));
    let room = Cube::new()
        .with_material(room_mat)
        .with_transform(Matrix::scaling(5., 5., 5.) * Matrix::translation(0., 1., 0.));
    objects.push(Rc::new(room));

    let table_top_pattern =
        Stripe::new(LBROWN, DBROWN).with_transform(Matrix::scaling(0.02, 0.02, 0.02));
    let table_top_mat = Material::default().with_pattern(table_top_pattern);
    let table_top = Cube::new()
        .with_material(table_top_mat)
        .with_transform(Matrix::translation(0., 1., 0.) * Matrix::scaling(1., 0.1, 1.));
    objects.push(Rc::new(table_top));

    let leg_mat_pattern =
        Stripe::new(LBROWN, DBROWN).with_transform(Matrix::scaling(0.02, 0.02, 0.02));
    let leg_mat = Material::default().with_pattern(leg_mat_pattern);
    let leg1 = Cube::new()
        .with_material(leg_mat)
        .with_transform(Matrix::translation(0.9, 0., 0.9) * Matrix::scaling(0.1, 1., 0.1));
    objects.push(Rc::new(leg1));

    let leg_mat_pattern =
        Stripe::new(LBROWN, DBROWN).with_transform(Matrix::scaling(0.02, 0.02, 0.02));
    let leg_mat = Material::default().with_pattern(leg_mat_pattern);
    let leg2 = Cube::new()
        .with_material(leg_mat)
        .with_transform(Matrix::translation(-0.9, 0., 0.9) * Matrix::scaling(0.1, 1., 0.1));
    objects.push(Rc::new(leg2));

    let leg_mat_pattern =
        Stripe::new(LBROWN, DBROWN).with_transform(Matrix::scaling(0.02, 0.02, 0.02));
    let leg_mat = Material::default().with_pattern(leg_mat_pattern);
    let leg3 = Cube::new()
        .with_material(leg_mat)
        .with_transform(Matrix::translation(-0.9, 0., -0.9) * Matrix::scaling(0.1, 1., 0.1));
    objects.push(Rc::new(leg3));

    let leg_mat_pattern =
        Stripe::new(LBROWN, DBROWN).with_transform(Matrix::scaling(0.02, 0.02, 0.02));
    let leg_mat = Material::default().with_pattern(leg_mat_pattern);
    let leg4 = Cube::new()
        .with_material(leg_mat)
        .with_transform(Matrix::translation(0.9, 0., -0.9) * Matrix::scaling(0.1, 1., 0.1));
    objects.push(Rc::new(leg4));

    let sphere_pattern = Solid::new(BLACK);
    let sphere_mat = Material::default()
        .with_diffuse(0.)
        .with_pattern(sphere_pattern)
        .with_transparency(1.)
        .with_refractive_index(1.5);
    let sphere = Sphere::new()
        .with_material(sphere_mat)
        .with_transform(Matrix::translation(0.5, 1.3, 0.) * Matrix::scaling(0.2, 0.2, 0.2));
    objects.push(Rc::new(sphere));

    let mirror_pattern = Solid::new(BLACK);
    let mirror_mat = Material::default()
        .with_diffuse(0.9)
        .with_pattern(mirror_pattern)
        .with_reflective(0.3);
    let mirror = Cube::new()
        .with_material(mirror_mat)
        .with_transform(Matrix::translation(0., 2., 4.9) * Matrix::scaling(2., 1., 0.01));
    objects.push(Rc::new(mirror));

    let cube_pattern = Solid::new(Color::new(0.1, 0.3, 0.9));
    let cube_mat = Material::default()
        .with_pattern(cube_pattern)
        .with_transparency(0.1);
    let cube = Cube::new().with_material(cube_mat).with_transform(
        Matrix::translation(-0.5, 1.25, -0.5)
            * Matrix::scaling(0.2, 0.2, 0.2)
            * Matrix::rotation(Axis::Y, -PI / 3.),
    );
    objects.push(Rc::new(cube));

    let camera_tf = camera::view_transform(
        Point::new(1.5, 1.5, -3.5),
        Point::new(0., 1., 0.),
        Vector::new(0., 1., 0.),
    );
    let camera = Camera::new(1920, 1080, PI / 3.).with_transform(camera_tf);

    let light = PointLight::new(Point::new(-4., 4., -4.), Color::new(1., 1., 1.));
    let world = World::new(objects, light);
    let canvas = world.render(camera, 5);
    println!("{}", canvas.to_ppm());
}
