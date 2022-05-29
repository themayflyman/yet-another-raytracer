use crate::aabb::{surrounding_box, AxisAlignedBoundingBox};
use crate::camera::degrees_to_radians;
use crate::material::{Isotropic, Material};
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::Vec3;
use std::sync::Arc;

use rand::Rng;

pub trait Hittable: Send + Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AxisAlignedBoundingBox>;

    fn pdf_value(&self, _origin: Vec3, _direction: Vec3) -> f64 {
        0.0
    }

    fn random(&self, _origin: Vec3) -> Vec3 {
        Vec3::new(1.0, 0.0, 0.0)
    }
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
    pub objects: Vec<Arc<dyn Hittable>>,
}

impl HittableList {
    pub fn new() -> Self {
        let objects: Vec<Arc<dyn Hittable>> = Vec::new();
        Self { objects }
    }

    pub fn add_object(&mut self, object: Arc<dyn Hittable>) {
        self.objects.push(object);
    }

    pub fn size(&self) -> usize {
        self.objects.len()
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut hit_rec: Option<HitRecord> = None;
        let mut closet_so_far: f64 = t_max;

        for sphere in self.objects.iter() {
            if let Some(temp_rec) = sphere.hit(ray, t_min, closet_so_far) {
                closet_so_far = temp_rec.t;
                hit_rec = Some(temp_rec);
            }
        }

        hit_rec
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AxisAlignedBoundingBox> {
        if self.objects.is_empty() {
            return None;
        }

        let mut temp_box: AxisAlignedBoundingBox = AxisAlignedBoundingBox::default();
        let mut first_box = true;

        for sphere in self.objects.iter() {
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

    fn pdf_value(&self, origin: Vec3, direction: Vec3) -> f64 {
        let weight = 1.0 / self.objects.len() as f64;

        return self
            .objects
            .iter()
            .map(|object| weight * object.pdf_value(origin, direction))
            .sum();
    }

    fn random(&self, origin: Vec3) -> Vec3 {
        if self.objects.is_empty() {
            return Vec3::new(1.0, 0.0, 0.0);
        } else {
            return self.objects[rand::thread_rng().gen_range::<usize>(0, self.objects.len() - 1)]
                .random(origin);
        }
    }
}

pub struct Translate<T: Hittable> {
    hittable: Arc<T>,
    offset: Vec3,
}

impl<T: Hittable> Translate<T> {
    pub fn new(hittable: Arc<T>, offset: Vec3) -> Self {
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

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AxisAlignedBoundingBox> {
        if let Some(output_box) = self.hittable.bounding_box(time0, time1) {
            return Some(AxisAlignedBoundingBox::new(
                output_box.min + self.offset,
                output_box.max + self.offset,
            ));
        }
        None
    }
}

pub struct RotateY<T: Hittable> {
    hittable: Arc<T>,
    sin_theta: f64,
    cos_theta: f64,
    bbox: Option<AxisAlignedBoundingBox>,
}

impl<T: Hittable> RotateY<T> {
    pub fn new(hittable: Arc<T>, angle: f64) -> Self {
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

                Some(AxisAlignedBoundingBox::new(min, max))
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

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AxisAlignedBoundingBox> {
        self.bbox
    }
}

pub struct ConstantMedium<TH: Hittable, TM: Material> {
    boundary: TH,
    phase_function: TM,
    neg_inv_density: f64,
}

impl<TT: Texture, TH: Hittable> ConstantMedium<TH, Isotropic<TT>> {
    pub fn new(boundary: TH, density: f64, texture: TT) -> Self {
        Self {
            boundary,
            phase_function: Isotropic::new(texture),
            neg_inv_density: -1.0 / density,
        }
    }
}

impl<TH: Hittable, TM: Material> Hittable for ConstantMedium<TH, TM> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut rng = rand::thread_rng();
        match self.boundary.hit(ray, -f64::INFINITY, f64::INFINITY) {
            Some(mut rec1) => match self.boundary.hit(ray, rec1.t + 0.0001, f64::INFINITY) {
                Some(mut rec2) => {
                    if rec1.t < t_min {
                        rec1.t = t_min;
                    }
                    if rec2.t > t_max {
                        rec2.t = t_max;
                    }
                    if rec1.t < rec2.t {
                        if rec1.t < 0.0 {
                            rec1.t = 0.0;
                        }

                        let ray_length = ray.direction().length();
                        let distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
                        let hit_distance = self.neg_inv_density * rng.gen::<f64>().ln();

                        if hit_distance < distance_inside_boundary {
                            let t = rec1.t + hit_distance / ray_length;
                            let p = ray.at(t);

                            Some(HitRecord {
                                u: 0.0,
                                v: 0.0,
                                t,
                                p,
                                normal: Vec3::new(1.0, 0.0, 0.0),
                                front_face: true,
                                material: &self.phase_function,
                            })
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                }

                None => None,
            },

            None => None,
        }
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AxisAlignedBoundingBox> {
        self.boundary.bounding_box(time0, time1)
    }
}

pub struct FlipFace<T: Hittable> {
    hittable: Arc<T>,
}

impl<T: Hittable> FlipFace<T> {
    pub fn new(hittable: Arc<T>) -> Self {
        Self { hittable }
    }
}

impl<T: Hittable> Hittable for FlipFace<T> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.hittable.hit(ray, t_min, t_max).map(|rec| HitRecord {
            u: rec.u,
            v: rec.v,
            t: rec.t,
            p: rec.p,
            normal: rec.normal,
            front_face: !rec.front_face,
            material: rec.material,
        })
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AxisAlignedBoundingBox> {
        self.hittable.bounding_box(time0, time1)
    }
}
