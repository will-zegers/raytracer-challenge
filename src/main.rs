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
    let floor = Plane::new(Matrix::scaling(0.25, 0.25, 0.25), floor_mat);
    objects.push(Rc::new(floor));

    let mut m_outer = Material::default();
    m_outer.pattern = Box::new(Solid::new(BLACK));
    m_outer.ambient = 0.5;
    m_outer.specular = 0.;
    m_outer.shininess = 0.;
    m_outer.transparency = 1.;
    m_outer.refractive_index = 1.5;
    let outer = Sphere::new()
        .set_transform(Matrix::translation(0., 1., 0.))
        .set_material(m_outer);
    objects.push(Rc::new(outer));

    let mut m_inner = Material::default();
    m_inner.pattern = Box::new(Solid::new(BLACK));
    m_inner.specular = 0.;
    m_inner.shininess = 0.;
    m_inner.transparency = 1.;
    m_inner.refractive_index = 1.;
    let inner = Sphere::new()
        .set_transform(Matrix::translation(0., 1., 0.) * Matrix::scaling(0.5, 0.5, 0.5))
        .set_material(m_inner);
    objects.push(Rc::new(inner));

    let camera_tf = camera::view_transform(
        Point::new(0., 1.5, -5.),
        Point::new(0., 1., 0.),
        Vector::new(0., 1., 0.),
    );
    let camera = Camera::new(640, 480, PI / 3.).set_transform(camera_tf);

    let world = World::new(objects, PointLight::default());
    let canvas = world.render(camera, 5);
    println!("{}", canvas.to_ppm());
}
