use rand::Rng;

use crate::aabb::AxisAlignedBoundingBox;
use crate::hittable::{HitRecord, Hittable};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct XYRect<M: Material> {
    pub x0: f32,
    pub x1: f32,
    pub y0: f32,
    pub y1: f32,
    pub k: f32,
    pub material: M,
}

impl<M: Material> XYRect<M> {
    pub fn new(x0: f32, x1: f32, y0: f32, y1: f32, k: f32, material: M) -> Self {
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

impl<M: Material> Hittable for XYRect<M> {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        let outputbox = AxisAlignedBoundingBox::new(
            Vec3::new(self.x0, self.y0, self.k - 0.0001),
            Vec3::new(self.x1, self.y1, self.k + 0.0001),
        );
        Some(outputbox)
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
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

        Some(HitRecord {
            u,
            v,
            t,
            p,
            normal,
            front_face,
            material: &self.material,
        })
    }
}

#[derive(Clone)]
pub struct XZRect<M: Material> {
    pub x0: f32,
    pub x1: f32,
    pub z0: f32,
    pub z1: f32,
    pub k: f32,
    pub material: M,
}

impl<M: Material> XZRect<M> {
    pub fn new(x0: f32, x1: f32, z0: f32, z1: f32, k: f32, material: M) -> Self {
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

impl<M: Material> Hittable for XZRect<M> {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        let outputbox = AxisAlignedBoundingBox::new(
            Vec3::new(self.x0, self.k - 0.0001, self.z0),
            Vec3::new(self.x1, self.k + 0.0001, self.z1),
        );
        Some(outputbox)
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
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

        Some(HitRecord {
            u,
            v,
            t,
            p,
            normal,
            front_face,
            material: &self.material,
        })
    }

    fn pdf_value(&self, origin: Vec3, direction: Vec3, wavelength: f32) -> f32 {
        if let Some(rec) = self.hit(&Ray::new(origin, direction, 0.0, wavelength), 0.001, f32::INFINITY) {
            let area = (self.x1 - self.x0) * (self.z1 - self.z0);
            let distance_squared = rec.t * rec.t * direction.length_squared();
            let cosine = direction.dot(&rec.normal).abs() / direction.length();

            distance_squared / (cosine * area)
        } else {
            0.0
        }
    }

    fn random(&self, origin: Vec3) -> Vec3 {
        let random_point = Vec3::new(
            rand::thread_rng().gen_range(self.x0..self.x1),
            self.k,
            rand::thread_rng().gen_range(self.z0..self.z1),
        );
        random_point - origin
    }
}

#[derive(Clone)]
pub struct YZRect<M: Material> {
    pub y0: f32,
    pub y1: f32,
    pub z0: f32,
    pub z1: f32,
    pub k: f32,
    pub material: M,
}

impl<M: Material> YZRect<M> {
    pub fn new(y0: f32, y1: f32, z0: f32, z1: f32, k: f32, material: M) -> Self {
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

impl<M: Material> Hittable for YZRect<M> {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        let outputbox = AxisAlignedBoundingBox::new(
            Vec3::new(self.k - 0.0001, self.y0, self.z0),
            Vec3::new(self.k + 0.0001, self.y1, self.z1),
        );
        Some(outputbox)
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
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

        Some(HitRecord {
            u,
            v,
            t,
            p,
            normal,
            front_face,
            material: &self.material,
        })
    }
}
