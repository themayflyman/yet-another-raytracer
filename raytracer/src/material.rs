use std::f32::consts::PI;
use std::sync::Arc;

use rand::{random, Rng};

use crate::hittable::HitRecord;
use crate::pdf::{CosinePDF, Pdf};
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct Scatter {
    // pub color: RGB,
    pub attenuation: f32,
    pub ray: Option<Ray>,
    pub is_specular: bool,
    pub pdf: Option<Arc<dyn Pdf>>,
}

pub trait Material: Send + Sync {
    fn scatter(&self, _ray_in: &Ray, _hit_record: &HitRecord) -> Option<Scatter> {
        None
    }
    // fn emitted(&self, _rec: &HitRecord, _u: f32, _v: f32, _p: Vec3) -> RGB {
    fn emitted(&self, _ray_in: &Ray, _hit_record: &HitRecord) -> f32 {
        return 0.0;
    }
    fn scatter_pdf(&self, _ray_in: &Ray, _hit: &HitRecord, _scattered: &Ray) -> f32 {
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
    fn scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<Scatter> {
        Some(Scatter {
            // color: self.albedo.value(ray_in.u, hit_record),
            attenuation: self.albedo.value(ray_in, hit_record),
            ray: None,
            is_specular: false,
            pdf: Some(Arc::new(CosinePDF::new(hit_record.normal))),
        })
    }
    fn scatter_pdf(&self, _ray_in: &Ray, hit: &HitRecord, scattered: &Ray) -> f32 {
        let cosine = hit.normal.dot(&scattered.direction().unit_vector());
        if cosine < 0.0 {
            0.0
        } else {
            cosine / PI
        }
    }
}

#[derive(Clone)]
pub struct Metal<T: Texture> {
    albedo: T,
    fuzz: f32,
}

impl<T: Texture> Metal<T> {
    pub fn new(albedo: T, fuzz: f32) -> Self {
        Metal { albedo, fuzz }
    }
}

fn reflect(v: Vec3, n: Vec3) -> Vec3 {
    v - 2.0 * v.dot(&n) * n
}

impl<T: Texture> Material for Metal<T> {
    fn scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<Scatter> {
        let reflected = reflect(ray_in.direction().unit_vector(), hit_record.normal);
        let scattered = Ray::new(
            hit_record.p,
            reflected + self.fuzz * random_in_unit_sphere(),
            ray_in.time(),
            ray_in.wavelength,
        );

        Some(Scatter {
            attenuation: self.albedo.value(ray_in, hit_record),
            ray: Some(scattered),
            is_specular: true,
            pdf: None,
        })
    }
}

// fn refract(uv: Vec3, n: Vec3, etai_over_etat: f32) -> Vec3 {
//     let cos_theta = f32::min(-uv.dot(&n), 1.0);
//     let r_out_perp = etai_over_etat * (uv + cos_theta * n);
//     let r_out_parallel = -(f32::abs(1.0 - r_out_perp.length_squared())).sqrt() * n;
//     r_out_perp + r_out_parallel
// }

fn reflectance(cosine: f32, ref_idx: f32) -> f32 {
    let mut r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}

#[derive(Clone, Copy)]
pub struct Dielectric {
    b1: f32,
    b2: f32,
    b3: f32,
    c1: f32,
    c2: f32,
    c3: f32,
}

// https://refractiveindex.info/?shelf=glass&book=BAF10&page=SCHOTT
pub static BAF10: Dielectric = Dielectric {
    b1: 1.5851495,
    b2: 0.143559385,
    b3: 1.08521269,
    c1: 0.00926681282 * 1e6,
    c2: 0.0424489805 * 1e6,
    c3: 105.613573 * 1e6,
};

// https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
pub static BK7: Dielectric = Dielectric {
    b1: 1.03961212,
    b2: 0.231792344,
    b3: 1.01046945,
    c1: 0.00600069867,
    c2: 0.0200179144,
    c3: 103.560653,
};

// https://refractiveindex.info/?shelf=glass&book=SF11&page=SCHOTT
pub static SF11: Dielectric = Dielectric {
    b1: 1.73759695,
    b2: 0.313747346,
    b3: 1.89878101,
    c1: 0.013188707,
    c2: 0.0623068142,
    c3: 155.23629,
};

// https://refractiveindex.info/?shelf=glass&book=FK51Apage=SCHOTT
pub static FK51A: Dielectric = Dielectric {
    b1: 0.971247817,
    b2: 0.216901417,
    b3: 0.904651666,
    c1: 0.00472301995,
    c2: 0.0153575612,
    c3: 168.68133,
};

// https://refractiveindex.info/?shelf=glass&book=LASF9&page=SCHOTT
pub static LASF9: Dielectric = Dielectric {
    b1: 2.00029547,
    b2: 0.298926886,
    b3: 1.80691843,
    c1: 0.0121426017,
    c2: 0.0538736236,
    c3: 156.530829,
};

// https://refractiveindex.info/?shelf=glass&book=SCHOTT-SF&page=N-SF66
pub static SF66: Dielectric = Dielectric {
    b1: 2.0245976,
    b2: 0.470187196,
    b3: 2.59970433,
    c1: 0.0147053225,
    c2: 0.0692998276,
    c3: 161.817601,
};

// impl Dielectric {
//     // pub fn new(index_of_refraction: f32) -> Self {
//     //     Dielectric {
//     //         index_of_refraction,
//     //     }
//     // }
// }

fn refract(v: Vec3, n: Vec3, ni_over_nt: f32) -> Option<Vec3> {
    let uv = v.unit_vector();
    let dt = uv.dot(&n);
    let discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0 - dt * dt);
    if discriminant > 0.0 {
        let refracted = (uv - n * dt) * ni_over_nt - n * f32::sqrt(discriminant);
        Some(refracted)
    } else {
        None
    }
}

fn schlick(cosine: f32, ref_idx: f32) -> f32 {
    let r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    let r0 = r0 * r0;
    r0 + (1.0 - r0) * f32::powi(1.0 - cosine, 5)
}

impl Material for Dielectric {
    fn scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<Scatter> {
        // let refraction_ratio = if hit_record.front_face {
        //     1.0 / self.index_of_refraction
        // } else {
        //     self.index_of_refraction
        // };

        // let unit_direction = ray_in.direction().unit_vector();
        // let cos_theta = f32::min(-ray_in.direction().unit_vector().dot(&hit_record.normal), 1.0);
        // let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        // let cannot_refract = refraction_ratio * sin_theta > 1.0;

        // if cannot_refract || reflectance(cos_theta, refraction_ratio) > random::<f32>() {
        //     Some(Scatter {
        //         attenuation: 1.0,
        //         ray: Some(Ray::new(
        //             hit_record.p,
        //             reflect(unit_direction, hit_record.normal),
        //             ray_in.time(),
        //             ray_in.wavelength,
        //         )),
        //         is_specular: true,
        //         pdf: None,
        //     })
        // } else {
        //     Some(Scatter {
        //         // color: Vec3::new(1.0, 1.0, 1.0),
        //         attenuation: 1.0,
        //         ray: Some(Ray::new(
        //             hit_record.p,
        //             refract(unit_direction, hit_record.normal, refraction_ratio),
        //             ray_in.time(),
        //             ray_in.wavelength,
        //         )),
        //         is_specular: true,
        //         pdf: None,
        //     })
        // }
        let wl_pow2 = ray_in.wavelength * ray_in.wavelength;
        // Sellmeier equation: https://en.wikipedia.org/wiki/Sellmeier_equation
        let refractive_index_squared = 1.0
            + self.b1 * wl_pow2 / (wl_pow2 - self.c1)
            + self.b2 * wl_pow2 / (wl_pow2 - self.c2)
            + self.b3 * wl_pow2 / (wl_pow2 - self.c3);
        let refractive_index = refractive_index_squared.sqrt();
        let (outward_normal, ni_over_nt, cosine) = if ray_in.direction().dot(&hit_record.normal)
            > 0.0
        {
            (
                -hit_record.normal,
                refractive_index,
                refractive_index * ray_in.direction().dot(&hit_record.normal) / ray_in.direction().length(),
            )
        } else {
            (
                hit_record.normal,
                refractive_index.recip(),
                -ray_in.direction().dot(&hit_record.normal) / ray_in.direction().length(),
            )
        };

        if let Some(refracted) = refract(ray_in.direction(), outward_normal, ni_over_nt) {
            let ray_out: Ray;
            if random::<f32>() < schlick(cosine, refractive_index) {
                let reflected = reflect(ray_in.direction(), hit_record.normal);
                ray_out = Ray::new(hit_record.p, reflected,  ray_in.time(), ray_in.wavelength);
            } else {
                ray_out = Ray::new(hit_record.p, refracted, ray_in.time(), ray_in.wavelength);
            }
            return Some(Scatter {
                attenuation: 1.0,
                is_specular: true,
                ray: Some(ray_out),
                pdf: None,
            });
        } else {
            let reflected = reflect(ray_in.direction(), hit_record.normal);
            return Some(Scatter {
                attenuation: 1.0,
                is_specular: true,
                ray: Some(Ray::new(
                    hit_record.p,
                    reflected,
                    ray_in.time(),
                    ray_in.wavelength,
                )),
                pdf: None,
            });
        }
    }
}

#[allow(dead_code)]
pub fn random_in_unit_vector() -> Vec3 {
    random_in_unit_sphere().unit_vector()
}

pub fn random_in_unit_sphere() -> Vec3 {
    const MAX: f32 = 1.0;
    const MIN: f32 = -1.0;
    let mut rng = rand::thread_rng();

    loop {
        let p: Vec3 = Vec3::new(
            rng.gen_range::<f32>(MIN, MAX),
            rng.gen_range::<f32>(MIN, MAX),
            rng.gen_range::<f32>(MIN, MAX),
        );
        if p.length_squared() >= 1.0 {
            continue;
        }
        return p;
    }
}
//
// #[allow(dead_code)]
// pub fn random_in_hemisphere(normal: &Vec3) -> Vec3 {
//     let in_unit_sphere = random_in_unit_sphere();
//     if in_unit_sphere.dot(normal) > 0.0 {
//         in_unit_sphere
//     } else {
//         -in_unit_sphere
//     }
// }
//
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
    fn emitted(&self, ray_in: &Ray, hit_record: &HitRecord) -> f32 {
        if hit_record.front_face {
            self.emit.value(ray_in, hit_record)
        } else {
            return 0.0;
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
    fn scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<Scatter> {
        Some(Scatter {
            attenuation: self.albedo.value(ray_in, hit_record),
            ray: Some(Ray::new(
                hit_record.p,
                random_in_unit_sphere(),
                ray_in.time(),
                ray_in.wavelength,
            )),
            is_specular: false,
            pdf: None,
        })
    }
}

#[derive(Clone)]
pub struct NoMaterial;

impl Material for NoMaterial {}
