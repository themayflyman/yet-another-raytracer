use crate::aabb::{surrounding_box, AABB};
use crate::hittable::{HitRecord, Hittable};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct XYRect {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub k: f64,
    pub material: Box<dyn Material>,
}

impl XYRect {
    pub fn new(x0: f64, x1: f64, y0: f64, y1: f64, k: f64, material: Box<dyn Material>) -> Self {
        Self {
            x0,
            x1,
            y0,
            y1,
            k,
            material,
        }
    }
}

impl Hittable for XYRect {
    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AABB> {
        let outputbox = AABB::new(
            Vec3::new(self.x0, self.y0, self.k - 0.0001),
            Vec3::new(self.x1, self.y1, self.k + 0.0001),
        );
        Some(outputbox)
    }

    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin().z()) / ray.direction().z();
        if t < t_min || t > t_max {
            return None;
        }

        let x = ray.origin().x() + t * ray.direction().x();
        let y = ray.origin().y() + t * ray.direction().y();
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }

        let u = (x - self.x0) / (self.x1 - self.x0);
        let v = (y - self.y0) / (self.y1 - self.y0);
        let p = ray.at(t);
        let normal: Vec3;
        let front_face: bool;
        let outward_normal = Vec3::new(0.0, 0.0, 1.0);
        if ray.direction().dot(&outward_normal) < 0.0 {
            normal = outward_normal;
            front_face = true;
        } else {
            normal = -outward_normal;
            front_face = false;
        }

        return Some(HitRecord::new(
            u,
            v,
            t,
            p,
            normal,
            front_face,
            self.material.clone(),
        ));
    }
}

#[derive(Clone)]
pub struct XZRect {
    pub x0: f64,
    pub x1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
    pub material: Box<dyn Material>,
}

impl XZRect {
    pub fn new(x0: f64, x1: f64, z0: f64, z1: f64, k: f64, material: Box<dyn Material>) -> Self {
        Self {
            x0,
            x1,
            z0,
            z1,
            k,
            material,
        }
    }
}

impl Hittable for XZRect {
    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AABB> {
        let outputbox = AABB::new(
            Vec3::new(self.x0, self.k - 0.0001, self.z0),
            Vec3::new(self.x1, self.k + 0.0001, self.z1),
        );
        Some(outputbox)
    }

    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin().y()) / ray.direction().y();
        if t < t_min || t > t_max {
            return None;
        }

        let x = ray.origin().x() + t * ray.direction().x();
        let z = ray.origin().z() + t * ray.direction().z();
        if x < self.x0 || x > self.x1 || z < self.z0 || z > self.z1 {
            return None;
        }

        let u = (x - self.x0) / (self.x1 - self.x0);
        let v = (z - self.z0) / (self.z1 - self.z0);
        let p = ray.at(t);
        let normal: Vec3;
        let front_face: bool;
        let outward_normal = Vec3::new(0.0, 1.0, 0.0);
        if ray.direction().dot(&outward_normal) < 0.0 {
            normal = outward_normal;
            front_face = true;
        } else {
            normal = -outward_normal;
            front_face = false;
        }

        return Some(HitRecord::new(
            u,
            v,
            t,
            p,
            normal,
            front_face,
            self.material.clone(),
        ));
    }
}

#[derive(Clone)]
pub struct YZRect {
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
    pub material: Box<dyn Material>,
}

impl YZRect {
    pub fn new(y0: f64, y1: f64, z0: f64, z1: f64, k: f64, material: Box<dyn Material>) -> Self {
        Self {
            y0,
            y1,
            z0,
            z1,
            k,
            material,
        }
    }
}

impl Hittable for YZRect {
    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AABB> {
        let outputbox = AABB::new(
            Vec3::new(self.k - 0.0001, self.y0, self.z0),
            Vec3::new(self.k + 0.0001, self.y1, self.z1),
        );
        Some(outputbox)
    }

    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin().x()) / ray.direction().x();
        if t < t_min || t > t_max {
            return None;
        }

        let y = ray.origin().y() + t * ray.direction().y();
        let z = ray.origin().z() + t * ray.direction().z();
        if y < self.y0 || y > self.y1 || z < self.z0 || z > self.z1 {
            return None;
        }

        let u = (y - self.y0) / (self.y1 - self.y0);
        let v = (z - self.z0) / (self.z1 - self.z0);
        let p = ray.at(t);
        let normal: Vec3;
        let front_face: bool;
        let outward_normal = Vec3::new(1.0, 0.0, 0.0);
        if ray.direction().dot(&outward_normal) < 0.0 {
            normal = outward_normal;
            front_face = true;
        } else {
            normal = -outward_normal;
            front_face = false;
        }

        return Some(HitRecord::new(
            u,
            v,
            t,
            p,
            normal,
            front_face,
            self.material.clone(),
        ));
    }
}
