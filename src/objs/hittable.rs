use crate::data::ray::Ray;
use crate::data::vec3::Vec3;
use super::sphere::Sphere;

pub trait Hittable {
    fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rc: &mut HitRecord) -> bool;
}

pub struct HitRecord {
    pub t: f64,
    pub p: Vec3,
    pub normal: Vec3,
    pub front_face: bool,
}

impl HitRecord {
    pub fn new() -> HitRecord {
        HitRecord {
            t: 0.0,
            p: Vec3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 0.0),
            front_face: false
        }
    }

    pub fn set_face_normal(& mut self, ray: &Ray, outward_normal: Vec3) {
        let front_face: bool = ray.direction().dot(&outward_normal) < 0.0;
        if front_face {
            self.normal = outward_normal;
        } else {
            self.normal = -outward_normal;
        }
    }
}

pub struct HittableList {
    spheres: Vec<Sphere>,
}

impl HittableList {
    pub fn new() -> HittableList {
        let sphere_list: Vec<Sphere> = Vec::new();
        HittableList {spheres: sphere_list}
    }

    pub fn add_sphere(&mut self, sphere: Sphere) {
        self.spheres.push(sphere);
    }

    pub fn size(&self) -> usize {
        self.spheres.len()
    }
}

impl Hittable for HittableList {
    fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rec: &mut HitRecord) -> bool {
        let mut temp_rec: HitRecord = HitRecord::new();
        let mut hit_anything: bool = false;
        let mut closet_so_far: f64 = t_max;

        for sphere in self.spheres.iter() {
            if sphere.intersect(ray, t_min, closet_so_far, &mut temp_rec) {
                hit_anything = true;
                closet_so_far = temp_rec.t;
                rec.t = temp_rec.t;
                rec.p = temp_rec.p;
                rec.normal = temp_rec.normal;
            }
        }
        hit_anything
    }
}
