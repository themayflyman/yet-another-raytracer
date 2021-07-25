use super::hittable::{ Hittable, HitRecord };
use crate::data::ray::Ray;
use crate::data::vec3::Vec3;

pub struct Sphere {
    center: Vec3,
    radius: f64,
}

impl Sphere {
    pub fn new(center: Vec3, radius: f64) -> Sphere {
        Sphere { center, radius}
    }

    pub fn center(&self) -> Vec3 {
        self.center
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }
}

impl Hittable for Sphere {
    fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rec: &mut HitRecord) -> bool {
        let oc: Vec3 = ray.origin() - self.center();
        let a: f64 = ray.direction().length_squared();
        let half_b: f64 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius() * self.radius();

        let discriminant: f64 = half_b * half_b - a * c;
        if discriminant < 0.0 { return false; }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f64 = (0.0 - half_b - discriminant.sqrt()) / a;
        if t < t_min || t_max < t {
            t = (0.0 - half_b + discriminant.sqrt()) / a;
            if t < t_min || t_max < t { return false; }
        }

        rec.t = t;
        rec.p = ray.at(t);
        let outward_normal: Vec3 = (rec.p - self.center()) / self.radius();
        rec.set_face_normal(ray, outward_normal);
        true
    }
}