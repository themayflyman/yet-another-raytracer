use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::material::Lambertian;
use crate::material::Material;

pub trait Hittable {
    fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rc: &mut HitRecord) -> bool;

    fn material(&self) -> Box<dyn Material>;
}

pub struct HitRecord {
    pub t: f64,
    pub p: Vec3,
    pub normal: Vec3,
    pub front_face: bool,
    pub material: Box<dyn Material>,
}

impl HitRecord {
    pub fn new() -> HitRecord {
        HitRecord {
            t: 0.0,
            p: Vec3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 0.0),
            front_face: false,
            material: Box::new(Lambertian::new(Vec3::new(0.0, 0.0, 0.0))),
        }
    }
}

pub struct HittableList {
    spheres: Vec<Box<dyn Hittable>>,
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

    pub fn intersect(&self, ray: &Ray, t_min: f64, t_max: f64, rec: &mut HitRecord) -> bool {
        let mut temp_rec: HitRecord = HitRecord::new();
        let mut hit_anything: bool = false;
        let mut closet_so_far: f64 = t_max;

        for sphere in self.spheres.iter() {
            if sphere.intersect(ray, t_min, closet_so_far, &mut temp_rec) {
                hit_anything = true;
                closet_so_far = temp_rec.t;
                rec.t = temp_rec.t;
                rec.p = temp_rec.p;
                rec.front_face = temp_rec.front_face;
                rec.normal = temp_rec.normal;
                rec.material = sphere.material();
            }
        }
        hit_anything
    }
}
