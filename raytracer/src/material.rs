use std::f64::consts::PI;
use std::sync::Arc;

use rand::{random, Rng};

use crate::hittable::HitRecord;
use crate::pdf::{CosinePDF, Pdf};
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct Scatter {
    pub color: Vec3,
    pub ray: Option<Ray>,
    pub is_specular: bool,
    pub pdf: Option<Arc<dyn Pdf>>,
}

pub trait Material: Send + Sync {
    fn scatter(&self, _ray_in: &Ray, _hit: &HitRecord) -> Option<Scatter> {
        None
    }
    fn emitted(&self, _rec: &HitRecord, _u: f64, _v: f64, _p: Vec3) -> Vec3 {
        Vec3::default()
    }
    fn scatter_pdf(&self, _ray_in: &Ray, _hit: &HitRecord, _scattered: &Ray) -> f64 {
        0.0
    }
}

#[derive(Clone)]
pub struct Lambertian<T: Texture> {
    albedo: T,
}

impl<T: Texture> Lambertian<T> {
    pub fn new(albedo: T) -> Self {
        Self { albedo }
    }
}

impl<T: Texture> Material for Lambertian<T> {
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Option<Scatter> {
        Some(Scatter {
            color: self.albedo.value(rec.u, rec.v, rec.p),
            ray: None,
            is_specular: false,
            pdf: Some(Arc::new(CosinePDF::new(rec.normal))),
        })
    }
    fn scatter_pdf(&self, _ray_in: &Ray, hit: &HitRecord, scattered: &Ray) -> f64 {
        let cosine = hit.normal.dot(&scattered.direction().unit_vector());
        if cosine < 0.0 {
            0.0
        } else {
            cosine / PI
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
        Metal { albedo, fuzz }
    }
}

fn reflect(v: Vec3, n: Vec3) -> Vec3 {
    v - 2.0 * v.dot(&n) * n
}

impl Material for Metal {
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Option<Scatter> {
        let reflected = reflect(_r_in.direction().unit_vector(), rec.normal);
        let scattered = Ray::new(
            rec.p,
            reflected + self.fuzz * random_in_unit_sphere(),
            _r_in.time(),
        );

        Some(Scatter {
            color: self.albedo,
            ray: Some(scattered),
            is_specular: true,
            pdf: None,
        })
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
        Dielectric {
            index_of_refraction,
        }
    }
}

impl Material for Dielectric {
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Option<Scatter> {
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
            Some(Scatter {
                color: Vec3::new(1.0, 1.0, 1.0),
                ray: Some(Ray::new(
                    rec.p,
                    reflect(unit_direction, rec.normal),
                    _r_in.time(),
                )),
                is_specular: true,
                pdf: None,
            })
        } else {
            Some(Scatter {
                color: Vec3::new(1.0, 1.0, 1.0),
                ray: Some(Ray::new(
                    rec.p,
                    refract(unit_direction, rec.normal, refraction_ratio),
                    _r_in.time(),
                )),
                is_specular: true,
                pdf: None,
            })
        }
    }
}

#[allow(dead_code)]
pub fn random_in_unit_vector() -> Vec3 {
    random_in_unit_sphere().unit_vector()
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

#[allow(dead_code)]
pub fn random_in_hemisphere(normal: &Vec3) -> Vec3 {
    let in_unit_sphere = random_in_unit_sphere();
    if in_unit_sphere.dot(normal) > 0.0 {
        in_unit_sphere
    } else {
        -in_unit_sphere
    }
}

#[derive(Clone)]
pub struct DiffuseLight<T: Texture> {
    emit: T,
}

impl<T: Texture> DiffuseLight<T> {
    pub fn new(emit: T) -> Self {
        Self { emit }
    }
}

impl<T: Texture> Material for DiffuseLight<T> {
    fn emitted(&self, rec: &HitRecord, u: f64, v: f64, p: Vec3) -> Vec3 {
        if rec.front_face {
            self.emit.value(u, v, p)
        } else {
            Vec3::default()
        }
    }
}

#[derive(Clone)]
pub struct Isotropic<T: Texture> {
    albedo: T,
}

impl<T: Texture> Isotropic<T> {
    pub fn new(albedo: T) -> Self {
        Self { albedo }
    }
}

impl<T: Texture> Material for Isotropic<T> {
    fn scatter(&self, _ray_in: &Ray, _hit: &HitRecord) -> Option<Scatter> {
        Some(Scatter {
            color: self.albedo.value(_hit.u, _hit.v, _hit.p),
            ray: Some(Ray::new(_hit.p, random_in_unit_sphere(), _ray_in.time())),
            is_specular: false,
            pdf: None,
        })
    }
}

#[derive(Clone)]
pub struct NoMaterial;

impl Material for NoMaterial {}
