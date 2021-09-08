extern crate rand;

use std::f64::consts::PI;

use rand::Rng;

use crate::data::ray::Ray;
use crate::data::vec3::Vec3;

pub struct Camera {
    lower_left_corner: Vec3,
    horizontal: Vec3,
    vertical: Vec3,
    origin: Vec3,
    // lens_radius: f64,
    // u: Vec3,
    // v: Vec3,
    // w: Vec3,
}

fn random_in_unit_disk() -> Vec3 {
    let mut rng = rand::thread_rng();
    let mut p: Vec3;
    while {
        p = 2.0 * Vec3::new(rng.gen::<f64>(), rng.gen::<f64>(), 0.0) - Vec3::new(1.0, 1.0, 0.0);
        p.dot(&p) >= 1.0
    } {}
    p
}

impl Camera {
    // pub fn new(lookfrom: Vec3, lookat: Vec3, vup: Vec3, vfov: f64, aspect: f64, aperture: f64, focus_dist: f64) -> Camera {
    pub fn new() -> Camera {
        // let lens_radius: f64 = aperture / 2.0;
        //
        // let theta: f64 = vfov * PI / 100.0;
        // let half_height: f64 = (theta / 2.0).tan();
        // let half_width: f64 = aspect * half_height;
        //
        // let origin: Vec3 = lookfrom;
        //
        // let w: Vec3 = (lookfrom - lookat).unit_vector();
        // let u: Vec3 = vup.cross(&w).unit_vector();
        // let v: Vec3 = w.cross(&u);
        //
        // let lower_left_corner: Vec3 = origin - half_height * focus_dist * u - half_height * focus_dist * v - focus_dist * w;
        // let horizontal: Vec3 = 2.0 * half_width * focus_dist * u;
        // let vertical: Vec3 = 2.0 * half_height * focus_dist * v;
        //
        // Camera {
        //     lower_left_corner,
        //     horizontal,
        //     vertical,
        //     origin,
        //     lens_radius,
        //     u,
        //     v,
        //     w,
        // }

        const ASPECT_RATIO: f64 = 16.0 / 9.0;
        const VIEWPORT_HEIGHT: f64 = 2.0;
        const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
        const FOCAL_LENGTH: f64 = 1.0;

        let origin: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let horizontal: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0.0, 0.0);
        let vertical: Vec3 = Vec3::new(0.0, VIEWPORT_HEIGHT, 0.0);
        let lower_left_corner: Vec3 =
            origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, FOCAL_LENGTH);

        Camera {
            lower_left_corner,
            horizontal,
            vertical,
            origin,
        }
    }

    pub fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin,
        )
    }

    // pub fn get_ray(&self, u: f64, v: f64) -> Ray {
    //     let rd: Vec3 = self.lens_radius * random_in_unit_disk();
    //     let offset: Vec3 = self.u * rd.x() + self.v() * rd.y();
    //     Ray::new(
    //         self.origin + offset,
    //         self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin - offset
    //     )
    // }
}
