use super::vec3::Vec3;

pub struct Onb {
    pub u: Vec3,
    pub v: Vec3,
    pub w: Vec3,
}

impl Onb {
    pub fn build_from_w(normal: Vec3) -> Self {
        let w = normal.unit_vector();
        let a = if w.x().abs() > 0.9 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let v = w.cross(&a).unit_vector();
        let u = w.cross(&v);

        Self { u, v, w }
    }

    pub fn local(&self, a: Vec3) -> Vec3 {
        a.x() * self.u + a.y() * self.v + a.z() * self.w
    }
}
