use crate::hittable::Hittable;
use rand::Rng;
use std::f64::consts::PI;
use std::sync::Arc;

use super::onb::Onb;

use super::vec3::Vec3;

pub trait Pdf {
    fn value(&self, direction: Vec3) -> f64;
    fn generate(&self) -> Vec3;
}

pub fn random_cosine_direction() -> Vec3 {
    let r1 = rand::random::<f64>();
    let r2 = rand::random::<f64>();
    let z = (1.0 - r2).sqrt();

    let phi = 2.0 * PI * r1;
    let x = phi.cos() * r2.sqrt();
    let y = phi.sin() * r2.sqrt();

    Vec3::new(x, y, z)
}

pub struct CosinePDF {
    pub uvw: Onb,
}

impl CosinePDF {
    pub fn new(w: Vec3) -> Self {
        Self {
            uvw: Onb::build_from_w(w),
        }
    }
}

impl Pdf for CosinePDF {
    fn value(&self, direction: Vec3) -> f64 {
        let cosine = direction.unit_vector().dot(&self.uvw.w);
        if cosine <= 0.0 {
            0.0
        } else {
            cosine / PI
        }
    }

    fn generate(&self) -> Vec3 {
        self.uvw.local(random_cosine_direction())
    }
}

pub struct HittablePDF {
    pub origin: Vec3,
    pub hittable: Arc<dyn Hittable>,
}

impl HittablePDF {
    pub fn new(hittable: Arc<dyn Hittable>, origin: Vec3) -> Self {
        Self { origin, hittable }
    }
}

impl Pdf for HittablePDF {
    fn value(&self, direction: Vec3) -> f64 {
        self.hittable.pdf_value(self.origin, direction)
    }

    fn generate(&self) -> Vec3 {
        self.hittable.random(self.origin)
    }
}

pub struct MixurePDF {
    pub p0: Arc<dyn Pdf>,
    pub p1: Arc<dyn Pdf>,
}

impl MixurePDF {
    pub fn new(p0: Arc<dyn Pdf>, p1: Arc<dyn Pdf>) -> Self {
        Self { p0, p1 }
    }
}

impl Pdf for MixurePDF {
    fn value(&self, direction: Vec3) -> f64 {
        0.5 * self.p0.value(direction) + 0.5 * self.p1.value(direction)
    }

    fn generate(&self) -> Vec3 {
        if rand::thread_rng().gen_range::<f64>(0.0, 1.0) < 0.5 {
            self.p0.generate()
        } else {
            self.p1.generate()
        }
    }
}
