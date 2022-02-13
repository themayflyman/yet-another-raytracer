use crate::aabb::{surrounding_box, AABB};
use crate::camera::degrees_to_radians;
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

pub trait Hittable: Send + Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB>;
}

pub struct HitRecord<'a> {
    pub u: f64,
    pub v: f64,
    pub t: f64,
    pub p: Vec3,
    pub normal: Vec3,
    pub front_face: bool,
    pub material: &'a dyn Material,
}

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

pub struct Translate<T: Hittable> {
    hittable: T,
    offset: Vec3,
}

impl<T: Hittable> Translate<T> {
    pub fn new(hittable: T, offset: Vec3) -> Self {
        Self { hittable, offset }
    }
}

impl<T: Hittable> Hittable for Translate<T> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut hit_rec: Option<HitRecord> = None;
        let moved_ray = Ray::new(ray.origin() - self.offset, ray.direction(), ray.time());

        if let Some(mut temp_rec) = self.hittable.hit(&moved_ray, t_min, t_max) {
            temp_rec.p = temp_rec.p + self.offset;
            hit_rec = Some(temp_rec);
        }

        hit_rec
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB> {
        if let Some(output_box) = self.hittable.bounding_box(time0, time1) {
            return Some(AABB::new(
                output_box.min + self.offset,
                output_box.max + self.offset,
            ));
        }
        return None;
    }
}

pub struct RotateY<T: Hittable> {
    hittable: T,
    sin_theta: f64,
    cos_theta: f64,
    bbox: Option<AABB>,
}

impl<T: Hittable> RotateY<T> {
    pub fn new(hittable: T, angle: f64) -> Self {
        let radians = degrees_to_radians(angle);
        let sin_theta = f64::sin(radians);
        let cos_theta = f64::cos(radians);
        let bbox = match hittable.bounding_box(0.0, 1.0) {
            Some(bbox) => {
                let mut min = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
                let mut max = Vec3::new(-f64::INFINITY, -f64::INFINITY, -f64::INFINITY);

                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..2 {
                            let x = i as f64 * bbox.max.x() + (1.0 - i as f64) * bbox.min.x();
                            let y = j as f64 * bbox.max.y() + (1.0 - j as f64) * bbox.min.y();
                            let z = k as f64 * bbox.max.z() + (1.0 - k as f64) * bbox.min.z();

                            let newx = cos_theta * x + sin_theta * z;
                            let newz = -sin_theta * x + cos_theta * z;

                            let tester = Vec3::new(newx, y, newz);

                            for c in 0..3 {
                                min[c] = f64::min(min[c], tester[c]);
                                max[c] = f64::max(max[c], tester[c]);
                            }
                        }
                    }
                }

                Some(AABB::new(min, max))
            }

            None => None,
        };

        Self {
            hittable,
            sin_theta,
            cos_theta,
            bbox,
        }
    }
}

impl<T: Hittable> Hittable for RotateY<T> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut origin = ray.origin();
        let mut direction = ray.direction();

        origin[0] = self.cos_theta * ray.origin()[0] - self.sin_theta * ray.origin()[2];
        origin[2] = self.sin_theta * ray.origin()[0] + self.cos_theta * ray.origin()[2];

        direction[0] = self.cos_theta * ray.direction()[0] - self.sin_theta * ray.direction()[2];
        direction[2] = self.sin_theta * ray.direction()[0] + self.cos_theta * ray.direction()[2];

        let rotated_ray = Ray::new(origin, direction, ray.time());

        match self.hittable.hit(&rotated_ray, t_min, t_max) {
            Some(hit_rec) => {
                let mut p = hit_rec.p;
                let mut normal = hit_rec.normal;

                p[0] = self.cos_theta * hit_rec.p[0] + self.sin_theta * hit_rec.p[2];
                p[2] = -self.sin_theta * hit_rec.p[0] + self.cos_theta * hit_rec.p[2];

                normal[0] = self.cos_theta * hit_rec.normal[0] + self.sin_theta * hit_rec.normal[2];
                normal[2] =
                    -self.sin_theta * hit_rec.normal[0] + self.cos_theta * hit_rec.normal[2];

                Some(HitRecord {
                    p,
                    normal,
                    ..hit_rec
                })
            }

            None => None,
        }
    }

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AABB> {
        self.bbox
    }
}
