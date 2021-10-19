use crate::aabb::{surrounding_box, AABB};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

pub trait HittableClone {
    fn clone_hittable<'a>(&self) -> Box<dyn Hittable>;
}

pub trait Hittable: HittableClone {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB>;
}

impl<T> HittableClone for T
where
    T: Hittable + Clone + 'static,
{
    fn clone_hittable(&self) -> Box<dyn Hittable> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Hittable> {
    fn clone(&self) -> Self {
        self.clone_hittable()
    }
}


pub struct HitRecord {
    pub u: f64,
    pub v: f64,
    pub t: f64,
    pub p: Vec3,
    pub normal: Vec3,
    pub front_face: bool,
    pub material: Box<dyn Material>,
}

impl HitRecord {
    pub fn new(
        u: f64,
        v: f64,
        t: f64,
        p: Vec3,
        normal: Vec3,
        front_face: bool,
        material: Box<dyn Material>,
    ) -> HitRecord {
        HitRecord {
            u,
            v,
            t,
            p,
            normal,
            front_face,
            material,
        }
    }
}

#[derive(Clone)]
pub struct HittableList {
    pub spheres: Vec<Box<dyn Hittable>>,
}

impl HittableList {
    pub fn new() -> HittableList {
        let sphere_list: Vec<Box<dyn Hittable>> = Vec::new();
        HittableList {
            spheres: sphere_list,
        }
    }

    pub fn add_sphere(&mut self, sphere: Box<dyn Hittable>) {
        self.spheres.push(sphere);
    }

    pub fn size(&self) -> usize {
        self.spheres.len()
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut hit_rec: Option<HitRecord> = None;
        let mut closet_so_far: f64 = t_max;

        for sphere in self.spheres.iter() {
            if let Some(temp_rec) = sphere.hit(ray, t_min, closet_so_far) {
                closet_so_far = temp_rec.t;
                hit_rec = Some(temp_rec);
            }
        }

        hit_rec
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB> {
        if self.spheres.is_empty() {
            return None;
        }

        let mut temp_box: AABB = AABB::default();
        let mut first_box = true;

        for sphere in self.spheres.iter() {
            if let Some(tot_box) = sphere.bounding_box(time0, time1) {
                if !first_box {
                    temp_box = surrounding_box(temp_box, tot_box);
                } else {
                    temp_box = tot_box;
                    first_box = false;
                }
            }
        }

        Some(temp_box)
    }
}
