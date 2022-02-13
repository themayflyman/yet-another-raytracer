use crate::vec3::Vec3;

#[derive(Clone, Copy, Debug)]
pub struct Ray {
    origin: Vec3,    // the ray origin
    direction: Vec3, // the ray direction
    time: f64,
}

impl Ray {
    pub fn new(origin: Vec3, direction: Vec3, time: f64) -> Ray {
        Ray {
            origin,
            direction,
            time,
        }
    }

    pub fn origin(&self) -> Vec3 {
        self.origin
    }

    pub fn direction(&self) -> Vec3 {
        self.direction
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    pub fn at(&self, t: f64) -> Vec3 {
        self.origin + t * self.direction
    }
}

impl PartialEq for Ray {
    fn eq(&self, other: &Ray) -> bool {
        self.origin == other.origin && self.direction == other.direction
    }
}