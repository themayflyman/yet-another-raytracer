use std::f64::consts::PI;

use crate::aabb::{surrounding_box, AxisAlignedBoundingBox};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

use super::hittable::{HitRecord, Hittable};

pub struct StillSphere<M: Material> {
    pub center: Vec3,
    pub radius: f64,
    pub material: M,
}

impl<M: Material> StillSphere<M> {
    pub fn new(center: Vec3, radius: f64, material: M) -> Self {
        Self {
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

impl<M: Material> Hittable for StillSphere<M> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc: Vec3 = ray.origin() - self.center();
        let a: f64 = ray.direction().length_squared();
        let half_b: f64 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius() * self.radius();

        let discriminant: f64 = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f64 = (0.0 - half_b - discriminant.sqrt()) / a;
        if t < t_min || t_max < t {
            t = (0.0 - half_b + discriminant.sqrt()) / a;
            if t < t_min || t_max < t {
                return None;
            }
        }

        let p: Vec3 = ray.at(t);
        let outward_normal: Vec3 = (p - self.center()) / self.radius();
        let normal: Vec3;
        let front_face: bool;
        if ray.direction().dot(&outward_normal) < 0.0 {
            normal = outward_normal;
            front_face = true;
        } else {
            normal = -outward_normal;
            front_face = false;
        }
        let (u, v) = get_sphere_uv(outward_normal);

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

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AxisAlignedBoundingBox> {
        Some(AxisAlignedBoundingBox::new(
            self.center - Vec3::new(self.radius, self.radius, self.radius),
            self.center + Vec3::new(self.radius, self.radius, self.radius),
        ))
    }
}

pub struct MovingSphere<M: Material> {
    center0: Vec3,
    center1: Vec3,
    time0: f64,
    time1: f64,
    pub radius: f64,
    pub material: M,
}

impl<M: Material> MovingSphere<M> {
    pub fn new(
        center0: Vec3,
        center1: Vec3,
        time0: f64,
        time1: f64,
        radius: f64,
        material: M,
    ) -> Self {
        Self {
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

impl<M: Material> Hittable for MovingSphere<M> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc: Vec3 = ray.origin() - self.center(ray.time());
        let a: f64 = ray.direction().length_squared();
        let half_b: f64 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius * self.radius;

        let discriminant: f64 = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f64 = (0.0 - half_b - discriminant.sqrt()) / a;
        if t < t_min || t_max < t {
            t = (0.0 - half_b + discriminant.sqrt()) / a;
            if t < t_min || t_max < t {
                return None;
            }
        }

        let p: Vec3 = ray.at(t);
        let outward_normal: Vec3 = (p - self.center(ray.time())) / self.radius;
        let normal: Vec3;
        let front_face: bool;
        if ray.direction().dot(&outward_normal) < 0.0 {
            normal = outward_normal;
            front_face = true;
        } else {
            normal = -outward_normal;
            front_face = false;
        }
        let (u, v) = get_sphere_uv(outward_normal);

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

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AxisAlignedBoundingBox> {
        let box0: AxisAlignedBoundingBox = AxisAlignedBoundingBox::new(
            self.center(time0) - Vec3::new(self.radius, self.radius, self.radius),
            self.center(time0) + Vec3::new(self.radius, self.radius, self.radius),
        );
        let box1: AxisAlignedBoundingBox = AxisAlignedBoundingBox::new(
            self.center(time1) - Vec3::new(self.radius, self.radius, self.radius),
            self.center(time1) + Vec3::new(self.radius, self.radius, self.radius),
        );
        Some(surrounding_box(box0, box1))
    }
}

pub fn get_sphere_uv(p: Vec3) -> (f64, f64) {
    // p: a given point on the sphere of radius one, centered at the origin.

    let theta = f64::acos(-p.y());
    let phi = f64::atan2(-p.z(), p.x()) + PI;

    (phi / (2.0 * PI), theta / PI)
}
