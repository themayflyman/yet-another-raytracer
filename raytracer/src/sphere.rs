use std::f32::consts::PI;

use crate::aabb::{surrounding_box, AxisAlignedBoundingBox};
use crate::material::Material;
use crate::onb::Onb;
use crate::ray::Ray;
use crate::vec3::Vec3;

use super::hittable::{HitRecord, Hittable};

fn random_to_sphere(radius: f32, distance_squared: f32) -> Vec3 {
    let r1 = rand::random::<f32>();
    let r2 = rand::random::<f32>();
    let z = 1.0 + r2 * ((1.0 - radius * radius / distance_squared).sqrt() - 1.0);

    let phi = 2.0 * PI * r1;
    let x = phi.cos() * (1.0 - z * z).sqrt();
    let y = phi.sin() * (1.0 - z * z).sqrt();

    Vec3::new(x, y, z)
}

pub struct StillSphere<M: Material> {
    pub center: Vec3,
    pub radius: f32,
    pub material: M,
}

impl<M: Material> StillSphere<M> {
    pub fn new(center: Vec3, radius: f32, material: M) -> Self {
        Self {
            center,
            radius,
            material,
        }
    }

    pub fn center(&self) -> Vec3 {
        self.center
    }

    pub fn radius(&self) -> f32 {
        self.radius
    }
}

impl<M: Material> Hittable for StillSphere<M> {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let oc: Vec3 = ray.origin() - self.center();
        let a: f32 = ray.direction().length_squared();
        let half_b: f32 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius() * self.radius();

        let discriminant: f32 = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f32 = (0.0 - half_b - discriminant.sqrt()) / a;
        if t < t_min || t_max < t {
            t = (0.0 - half_b + discriminant.sqrt()) / a;
            if t < t_min || t_max < t {
                return None;
            }
        }

        let p: Vec3 = ray.at(t);
        let outward_normal: Vec3 = (p - self.center()) / self.radius().abs();
        let (normal, front_face) = if self.radius() < 0.0 {
            (-outward_normal, ray.direction().dot(&outward_normal) > 0.0)
        } else {
            (outward_normal, ray.direction().dot(&outward_normal) < 0.0)
        };
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

    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        Some(AxisAlignedBoundingBox::new(
            self.center - Vec3::new(self.radius, self.radius, self.radius),
            self.center + Vec3::new(self.radius, self.radius, self.radius),
        ))
    }

    fn pdf_value(&self, origin: Vec3, direction: Vec3, wavelength: f32) -> f32 {
        if let Some(_rec) = self.hit(&Ray::new(origin, direction, 0.0, wavelength), 0.001, f32::INFINITY) {
            let cos_theta_max = (1.0
                - self.radius() * self.radius() / (self.center() - origin).length_squared())
            .sqrt();
            let solid_angle = 2.0 * PI * (1.0 - cos_theta_max);

            1.0 / solid_angle
        } else {
            0.0
        }
    }

    fn random(&self, origin: Vec3) -> Vec3 {
        let direction = self.center() - origin;
        let distance_squared = direction.length_squared();
        let uvw = Onb::build_from_w(direction);

        uvw.local(random_to_sphere(self.radius(), distance_squared))
    }
}

pub struct MovingSphere<M: Material> {
    center0: Vec3,
    center1: Vec3,
    time0: f32,
    time1: f32,
    pub radius: f32,
    pub material: M,
}

impl<M: Material> MovingSphere<M> {
    pub fn new(
        center0: Vec3,
        center1: Vec3,
        time0: f32,
        time1: f32,
        radius: f32,
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

    pub fn center(&self, time: f32) -> Vec3 {
        self.center0
            + ((time - self.time0) / (self.time1 - self.time0)) * (self.center1 - self.center0)
    }
}

impl<M: Material> Hittable for MovingSphere<M> {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let oc: Vec3 = ray.origin() - self.center(ray.time());
        let a: f32 = ray.direction().length_squared();
        let half_b: f32 = oc.dot(&ray.direction());
        let c = oc.length_squared() - self.radius * self.radius;

        let discriminant: f32 = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        // Find the nearest root that lies in the acceptable range.
        let mut t: f32 = (0.0 - half_b - discriminant.sqrt()) / a;
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

    fn bounding_box(&self, time0: f32, time1: f32) -> Option<AxisAlignedBoundingBox> {
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

pub fn get_sphere_uv(p: Vec3) -> (f32, f32) {
    // p: a given point on the sphere of radius one, centered at the origin.

    let theta = f32::acos(-p.y());
    let phi = f32::atan2(-p.z(), p.x()) + PI;

    (phi / (2.0 * PI), theta / PI)
}
