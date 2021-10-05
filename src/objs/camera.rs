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

fn degrees_to_radians(degrees: f64) -> f64 {
    return degrees * PI / 180.0;
}

impl Camera {
    pub fn new(lookfrom: Vec3, lookat: Vec3, vup: Vec3, vfov: f64, aspect_ratio: f64) -> Camera {
        const FOCAL_LENGTH: f64 = 1.0;
        let theta: f64 = degrees_to_radians(vfov);
        let h: f64 = f64::tan(theta / 2.0);
        let viewport_height: f64 = 2.0 * h;
        let viewport_width: f64 = aspect_ratio * viewport_height;

        let w: Vec3 = (lookfrom - lookat).unit_vector();
        let u: Vec3 = vup.cross(&w).unit_vector();
        let v: Vec3 = w.cross(&u);

        let origin: Vec3 = lookfrom;
        let horizontal: Vec3 = viewport_width * u;
        let vertical: Vec3 = viewport_height * v;
        let lower_left_corner: Vec3 = origin - horizontal / 2.0 - vertical / 2.0 - w;

        Camera {
            lower_left_corner,
            horizontal,
            vertical,
            origin,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left_corner + s * self.horizontal + t * self.vertical - self.origin,
        )
    }
}
