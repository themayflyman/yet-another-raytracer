use std::sync::Arc;

use crate::camera::Camera;
use crate::color::{gen_wavelength, HasReflectance, RGB, XYZ};
use crate::ray::Ray;
// use crate::sampler::Sampler;
use crate::hittable::Hittable;
use crate::pdf::{HittablePDF, MixurePDF, Pdf};
use crate::scene::Scene;

pub trait Integrator: Send + Sync {
    // fn render(&self, scene: &Scene);
    fn li(&self, r: &Ray, scene: &Scene) -> XYZ;
    // fn integrate(
    //   ray_in: &Ray,
    //   world: &HittableList,
    //   lights: Arc<HittableList>,
    //   background_color: &RGB,
    //   depth: usize,
    // )
}

pub struct SamplerIntegrator {
    pub camera: Camera,
    // pub sampler: Sampler,
    max_bounces: usize,
}

impl Integrator for SamplerIntegrator {
    fn li(&self, ray_in: &Ray, scene: &Scene) -> XYZ {
        let mut throughput = 1.0;
        let mut ray_tracing = ray_in.clone();
        let mut bounces = 0;

        while bounces < self.max_bounces {
            if let Some(hit_record) = scene.aggregate.hit(&ray_tracing, 0.001, f64::INFINITY) {
                let emitted = hit_record.material.emitted(&ray_tracing, &hit_record);

                if let Some(scattered) = hit_record.material.scatter(&ray_tracing, &hit_record) {
                    if let Some(scattered_ray) = scattered.ray {
                        throughput *= scattered.attenuation;
                        ray_tracing = scattered_ray;
                    } else {
                        if scattered.pdf.is_none() {
                            panic!("Pdf not provided");
                        }

                        let mixture_pdf = if scene.lights.objects.is_empty() {
                            let light_pdf = scattered.pdf.unwrap();
                            MixurePDF::new(light_pdf.clone(), light_pdf)
                        } else {
                            let light_pdf = HittablePDF::new(scene.lights.clone(), hit_record.p);
                            MixurePDF::new(Arc::new(light_pdf), scattered.pdf.unwrap())
                        };
                        let wl = gen_wavelength(
                            ray_tracing.wavelength - 10.0,
                            ray_tracing.wavelength + 10.0,
                        );
                        let s =
                            Ray::new(hit_record.p, mixture_pdf.generate(), ray_tracing.time(), wl);
                        let pdf_val = mixture_pdf.value(s.direction(), wl);

                        throughput *= scattered.attenuation
                            * hit_record
                                .material
                                .scatter_pdf(&ray_tracing, &hit_record, &s)
                            / pdf_val;

                        ray_tracing = s;
                    }
                    bounces += 1;
                } else {
                    return XYZ::from_wavelength(ray_tracing.wavelength) * (emitted + throughput);
                }
            } else {
                return XYZ::from_wavelength(ray_tracing.wavelength)
                    * scene.background_color.reflect(ray_tracing.wavelength)
                    * throughput;
            }
        }

        return XYZ::from_wavelength(ray_tracing.wavelength) * throughput;
    }
}
