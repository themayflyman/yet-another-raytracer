use super::vec3::Vec3;

#[derive(Debug)]
pub struct Ray {
    a: Vec3, // the ray origin
    b: Vec3, // the ray direction
}

impl Ray {
    pub fn new(a: Vec3, b: Vec3) -> Ray {
        Ray { a, b }
    }

    pub fn origin(&self) -> Vec3 {
        self.a
    }

    pub fn direction(&self) -> Vec3 {
        self.b
    }

    pub fn at(&self, t: f64) -> Vec3 {
        self.a + t * self.b
    }
}

impl PartialEq for Ray {
    fn eq(&self, other: &Ray) -> bool {
        self.origin() == other.origin() && self.direction() == other.direction()
    }
}
