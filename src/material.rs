use rand::{random, Rng};

use crate::data::ray::Ray;
use crate::data::vec3::Vec3;
use crate::objs::hittable::HitRecord;

#[derive(Clone, Debug)]
pub struct Scatter {
    pub color: Vec3,
    pub ray: Option<Ray>,
}

// A workaround to implement Clone for Box
pub trait MaterialClone {
    fn clone_material<'a>(&self) -> Box<dyn Material>;
}

pub trait Material: MaterialClone {
    fn scatter(&self, ray_in: &Ray, hit: &HitRecord) -> Scatter;
}

impl<T> MaterialClone for T
where
    T: Material + Clone + 'static,
{
    fn clone_material(&self) -> Box<dyn Material> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Material> {
    fn clone(&self) -> Self {
        self.clone_material()
    }
}

#[derive(Clone)]
pub struct Lambertian {
    albedo: Vec3,
}

impl Lambertian {
    pub fn new(albedo: Vec3) -> Lambertian {
        return Lambertian { albedo };
    }
}

impl Material for Lambertian {
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Scatter {
        let target = rec.p + rec.normal + random_in_unit_sphere();
        Scatter {
            color: self.albedo,
            ray: Some(Ray::new(rec.p, target - rec.p)),
        }
    }
}

#[derive(Clone)]
pub struct Metal {
    albedo: Vec3,
    fuzz: f64,
}

impl Metal {
    pub fn new(albedo: Vec3, fuzz: f64) -> Metal {
        return Metal { albedo, fuzz };
    }
}

fn reflect(v: Vec3, n: Vec3) -> Vec3 {
    v - 2.0 * v.dot(&n) * n
}

impl Material for Metal {
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Scatter {
        let reflected = reflect(_r_in.direction().unit_vector(), rec.normal);
        let scattered = Ray::new(rec.p, reflected + self.fuzz * random_in_unit_sphere());

        Scatter {
            color: self.albedo,
            ray: if scattered.direction().dot(&rec.normal) <= 0.0 {
                None
            } else {
                Some(scattered)
            },
        }
    }
}

fn refract(uv: Vec3, n: Vec3, etai_over_etat: f64) -> Vec3 {
    let cos_theta = f64::min(-uv.dot(&n), 1.0);
    let r_out_perp = etai_over_etat * (uv + cos_theta * n);
    let r_out_parallel = -(f64::abs(1.0 - r_out_perp.length_squared())).sqrt() * n;
    r_out_perp + r_out_parallel
}

fn reflectance(cosine: f64, ref_idx: f64) -> f64 {
    let mut r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}

#[derive(Clone)]
pub struct Dielectric {
    pub index_of_refraction: f64,
}

impl Dielectric {
    pub fn new(index_of_refraction: f64) -> Dielectric {
        return Dielectric {
            index_of_refraction,
        };
    }
}

impl Material for Dielectric {
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Scatter {
        let refraction_ratio = if rec.front_face {
            1.0 / self.index_of_refraction
        } else {
            self.index_of_refraction
        };

        let unit_direction = _r_in.direction().unit_vector();
        let cos_theta = f64::min(-_r_in.direction().unit_vector().dot(&rec.normal), 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        let cannot_refract = refraction_ratio * sin_theta > 1.0;

        if cannot_refract || reflectance(cos_theta, refraction_ratio) > random::<f64>() {
            return Scatter {
                color: Vec3::new(1.0, 1.0, 1.0),
                ray: Some(Ray::new(rec.p, reflect(unit_direction, rec.normal))),
            };
        } else {
            return Scatter {
                color: Vec3::new(1.0, 1.0, 1.0),
                ray: Some(Ray::new(
                    rec.p,
                    refract(unit_direction, rec.normal, refraction_ratio),
                )),
            };
        }
    }
}

pub fn random_in_unit_sphere() -> Vec3 {
    const MAX: f64 = 1.0;
    const MIN: f64 = -1.0;
    let mut rng = rand::thread_rng();

    loop {
        let p: Vec3 = Vec3::new(
            rng.gen_range::<f64>(MIN, MAX),
            rng.gen_range::<f64>(MIN, MAX),
            rng.gen_range::<f64>(MIN, MAX),
        );
        if p.length_squared() >= 1.0 {
            continue;
        }
        return p;
    }
}