use crate::vec3::Vec3;

#[derive(Clone, Copy, Debug)]
pub struct Ray {
    origin: Vec3,    // the ray origin
    direction: Vec3, // the ray direction
    time: f32,
    pub wavelength: f32
}

impl Ray {
    pub fn new(origin: Vec3, direction: Vec3, time: f32, wavelength: f32) -> Ray {
        Ray {
            origin,
            direction,
            time,
            wavelength
        }
    }

    pub fn origin(&self) -> Vec3 {
        self.origin
    }

    pub fn direction(&self) -> Vec3 {
        self.direction
    }

    pub fn time(&self) -> f32 {
        self.time
    }

    pub fn at(&self, t: f32) -> Vec3 {
        self.origin + t * self.direction
    }
}

impl PartialEq for Ray {
    fn eq(&self, other: &Ray) -> bool {
        self.origin == other.origin && self.direction == other.direction
    }
}
