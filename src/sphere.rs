use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::material::Material;

use super::hittable::{HitRecord, Hittable};

pub struct StillSphere {
    center: Vec3,
    radius: f64,
    pub material: Box<dyn Material>,
}

impl StillSphere {
    pub fn new(center: Vec3, radius: f64, material: Box<dyn Material>) -> StillSphere {
        StillSphere {
            center,
            radius,
            material,
        }
    }

    pub fn center(&self) -> Vec3 {
        self.center
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }
}

impl Hittable for StillSphere {
    fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rec: &mut HitRecord) -> bool {
        let oc: Vec3 = ray.origin() - self.center();
        let a: f64 = ray.direction().length_squared();
        let half_b: f64 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius() * self.radius();

        let discriminant: f64 = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return false;
        }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f64 = (0.0 - half_b - discriminant.sqrt()) / a;
        if t < t_min || t_max < t {
            t = (0.0 - half_b + discriminant.sqrt()) / a;
            if t < t_min || t_max < t {
                return false;
            }
        }

        rec.t = t;
        rec.p = ray.at(t);
        rec.material = self.material.clone();
        let outward_normal: Vec3 = (rec.p - self.center()) / self.radius();
        if ray.direction().dot(&outward_normal) < 0.0 {
            rec.normal = outward_normal;
            rec.front_face = true;
        } else {
            rec.normal = -outward_normal;
            rec.front_face = false;
        }
        true
    }

    fn material(&self) -> Box<dyn Material> {
        self.material.clone()
    }
}

pub struct MovingSphere {
    center0: Vec3,
    center1: Vec3,
    time0: f64,
    time1: f64,
    pub radius: f64,
    pub material: Box<dyn Material>,
}

impl MovingSphere {
    pub fn new(
        center0: Vec3,
        center1: Vec3,
        time0: f64,
        time1: f64,
        radius: f64,
        material: Box<dyn Material>,
    ) -> MovingSphere {
        MovingSphere {
            center0,
            center1,
            time0,
            time1,
            radius,
            material,
        }
    }

    pub fn center(&self, time: f64) -> Vec3 {
        self.center0
            + ((time - self.time0) / (self.time1 - self.time0)) * (self.center1 - self.center0)
    }
}

impl Hittable for MovingSphere {
    fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rec: &mut HitRecord) -> bool {
        let oc: Vec3 = ray.origin() - self.center(ray.time());
        let a: f64 = ray.direction().length_squared();
        let half_b: f64 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius * self.radius;

        let discriminant: f64 = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return false;
        }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f64 = (0.0 - half_b - discriminant.sqrt()) / a;
        if t < t_min || t_max < t {
            t = (0.0 - half_b + discriminant.sqrt()) / a;
            if t < t_min || t_max < t {
                return false;
            }
        }

        rec.t = t;
        rec.p = ray.at(t);
        rec.material = self.material.clone();
        let outward_normal: Vec3 = (rec.p - self.center(ray.time())) / self.radius;
        if ray.direction().dot(&outward_normal) < 0.0 {
            rec.normal = outward_normal;
            rec.front_face = true;
        } else {
            rec.normal = -outward_normal;
            rec.front_face = false;
        }
        true
    }

    fn material(&self) -> Box<dyn Material> {
        self.material.clone()
    }
}
