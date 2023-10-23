use crate::color::XYZ;

pub trait Integrator: Send + Sync {
    fn integrate(
      ray_in: &Ray,
      world: &HittableList,
      lights: Arc<HittableList>,
      background_color: &RGB,
      depth: usize,
    )
}

pub struct MonteCarloIntegrator {
    max_bounces: usize
}

impl Integrator for MonteCarloIntegrator {
    fn integrate(&self,
          ray_in: &Ray,
          world: &HittableList,
          lights: Arc<HittableList>,
          background_color: &RGB,
          depth: usize,
        ) {
        let bounce = 0;
        let color: XYZ = XYZ::from_wavelength(ray_in.wavelength);
        let reflectance = 1.0;

        let mut spawned_ray = ray_in.clone();

        loop {
            if bounce >= self.max_bounces {
                break;
            }

           if let Some(hit_record) = world.hit(spawned_ray, 0.001, f32::INFINITY) {
               let emitted = hit_record.material.emitted(ray_in, &hit_record);

               if let Some(scattered) = hit_record.material.scatter(ray_in, &hit_record) {
                   if let Some(scattered_ray) = scattered.ray {
                       reflectance += emitted
                       return emitted
                           + scattered.attenuation
                               * ray_reflectance(
                                   &scattered_ray,
                                   world,
                                   lights,
                                   background_color,
                                   depth - 1,
                               );
                   }
                   if scattered.pdf.is_none() {
                       panic!("Pdf not provided")
                   }

                   let mixure_pdf = if lights.objects.is_empty() {
                       let light_pdf = scattered.pdf.unwrap();
                       MixurePDF::new(light_pdf.clone(), light_pdf)
                   } else {
                       let light_pdf = HittablePDF::new(lights.clone(), hit_record.p);
                       MixurePDF::new(Arc::new(light_pdf), scattered.pdf.unwrap())
                   };
                   let wl = gen_wavelength(ray_in.wavelength - 10.0, ray_in.wavelength + 10.0);
                   let s = Ray::new(hit_record.p, mixure_pdf.generate(), ray_in.time(), wl);
                   let pdf_val = mixure_pdf.value(s.direction(), wl);
                   return emitted
                       + scattered.attenuation
                           * ray_reflectance(&s, world, lights, background_color, depth - 1)
                           * hit_record.material.scatter_pdf(ray_in, &hit_record, &s)
                           / pdf_val;
               } else {
                   return emitted;
               }
           }

           return background_color.reflect(ray_in.wavelength);
        }



    }
}

fn ray_color(
    r: &Ray,
    world: &HittableList,
    lights: Arc<HittableList>,
    background_color: &RGB,
    max_depth: usize,
) -> XYZ {
    let reflectance = ray_reflectance(r, world, lights, background_color, max_depth);
    return XYZ::from_wavelength(r.wavelength) * reflectance;
}

fn ray_reflectance(
    ray_in: &Ray,
    world: &HittableList,
    lights: Arc<HittableList>,
    background_color: &RGB,
    depth: usize,
) -> f32 {
    if depth == 0 {
        return 1.0;
    }

    if let Some(hit_record) = world.hit(ray_in, 0.001, f32::INFINITY) {
        let emitted = hit_record.material.emitted(ray_in, &hit_record);

        if let Some(scattered) = hit_record.material.scatter(ray_in, &hit_record) {
            if let Some(scattered_ray) = scattered.ray {
                return emitted
                    + scattered.attenuation
                        * ray_reflectance(
                            &scattered_ray,
                            world,
                            lights,
                            background_color,
                            depth - 1,
                        );
            }
            if scattered.pdf.is_none() {
                panic!("Pdf not provided")
            }

            let mixure_pdf = if lights.objects.is_empty() {
                let light_pdf = scattered.pdf.unwrap();
                MixurePDF::new(light_pdf.clone(), light_pdf)
            } else {
                let light_pdf = HittablePDF::new(lights.clone(), hit_record.p);
                MixurePDF::new(Arc::new(light_pdf), scattered.pdf.unwrap())
            };
            let wl = gen_wavelength(ray_in.wavelength - 10.0, ray_in.wavelength + 10.0);
            let s = Ray::new(hit_record.p, mixure_pdf.generate(), ray_in.time(), wl);
            let pdf_val = mixure_pdf.value(s.direction(), wl);
            return emitted
                + scattered.attenuation
                    * ray_reflectance(&s, world, lights, background_color, depth - 1)
                    * hit_record.material.scatter_pdf(ray_in, &hit_record, &s)
                    / pdf_val;
        } else {
            return emitted;
        }
    }

    return background_color.reflect(ray_in.wavelength);
}
