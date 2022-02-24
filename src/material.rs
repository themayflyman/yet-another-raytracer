use rand::{random, Rng};

use crate::hittable::HitRecord;
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::Vec3;

#[derive(Clone, Debug)]
pub struct Scatter {
    pub color: Vec3,
    pub ray: Option<Ray>,
}

pub trait Material: Send + Sync {
    fn scatter(&self, _ray_in: &Ray, _hit: &HitRecord) -> Scatter {
        Scatter {
            color: Vec3::default(),
            ray: None,
        }
    }
    fn emitted(&self, _u: f64, _v: f64, _p: Vec3) -> Vec3 {
        Vec3::default()
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
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Scatter {
        let mut scatter_direction = rec.normal + random_in_unit_vector();

        if scatter_direction.near_zero() {
            scatter_direction = rec.normal;
        }

        Scatter {
            color: self.albedo.value(rec.u, rec.v, rec.p),
            ray: Some(Ray::new(rec.p, scatter_direction, _r_in.time())),
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
    fn scatter(&self, _r_in: &Ray, rec: &HitRecord) -> Scatter {
        let reflected = reflect(_r_in.direction().unit_vector(), rec.normal);
        let scattered = Ray::new(
            rec.p,
            reflected + self.fuzz * random_in_unit_sphere(),
            _r_in.time(),
        );

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
        Dielectric {
            index_of_refraction,
        }
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
            Scatter {
                color: Vec3::new(1.0, 1.0, 1.0),
                ray: Some(Ray::new(
                    rec.p,
                    reflect(unit_direction, rec.normal),
                    _r_in.time(),
                )),
            }
        } else {
            Scatter {
                color: Vec3::new(1.0, 1.0, 1.0),
                ray: Some(Ray::new(
                    rec.p,
                    refract(unit_direction, rec.normal, refraction_ratio),
                    _r_in.time(),
                )),
            }
        }
    }
}

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
    fn scatter(&self, _ray_in: &Ray, _hit: &HitRecord) -> Scatter {
        Scatter {
            color: Vec3::default(),
            ray: None,
        }
    }

    fn emitted(&self, u: f64, v: f64, p: Vec3) -> Vec3 {
        self.emit.value(u, v, p)
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
    fn scatter(&self, _ray_in: &Ray, _hit: &HitRecord) -> Scatter {
        Scatter {
            color: self.albedo.value(_hit.u, _hit.v, _hit.p),
            ray: Some(Ray::new(_hit.p, random_in_unit_sphere(), _ray_in.time())),
        }
    }
}
