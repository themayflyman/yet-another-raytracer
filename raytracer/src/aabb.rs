use crate::ray::Ray;
use crate::vec3::Vec3;
use std::mem::swap;

#[derive(Default, Clone, Copy)]
pub struct AxisAlignedBoundingBox {
    pub min: Vec3,
    pub max: Vec3,
}

impl AxisAlignedBoundingBox {
    pub fn new(min: Vec3, max: Vec3) -> Self {
        Self { min, max }
    }

    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        let mut t_min = t_min;
        let mut t_max = t_max;

        for a in 0..3 {
            let inv_d: f32 = 1.0 / ray.direction()[a];
            let mut t0: f32 = (self.min[a] - ray.origin()[a]) * inv_d;
            let mut t1: f32 = (self.max[a] - ray.origin()[a]) * inv_d;
            if inv_d < 0.0 {
                swap(&mut t0, &mut t1);
            }
            t_min = f32::min(t0, t_min);
            t_max = f32::max(t1, t_max);
            // t_min = if t0 > t_min { t0 } else { t_min };
            // t_max = if t1 < t_max { t1 } else { t_max };
            if t_max <= t_min {
                return false;
            }
        }
        true
    }

    pub fn surrounding_box(box0: &AxisAlignedBoundingBox, box1: &AxisAlignedBoundingBox) -> Self {
        let small: Vec3 = Vec3::new(
            f32::min(box0.min.x(), box1.min.x()),
            f32::min(box0.min.y(), box1.min.y()),
            f32::min(box0.min.z(), box1.min.z()),
        );
        let big: Vec3 = Vec3::new(
            f32::max(box0.max.x(), box1.max.x()),
            f32::max(box0.max.y(), box1.max.y()),
            f32::max(box0.max.z(), box1.max.z()),
        );
        Self::new(small, big)
    }
}

pub fn surrounding_box(
    box0: AxisAlignedBoundingBox,
    box1: AxisAlignedBoundingBox,
) -> AxisAlignedBoundingBox {
    let small: Vec3 = Vec3::new(
        f32::min(box0.min.x(), box1.min.x()),
        f32::min(box0.min.y(), box1.min.y()),
        f32::min(box0.min.z(), box1.min.z()),
    );
    let big: Vec3 = Vec3::new(
        f32::max(box0.max.x(), box1.max.x()),
        f32::max(box0.max.y(), box1.max.y()),
        f32::max(box0.max.z(), box1.max.z()),
    );
    AxisAlignedBoundingBox::new(small, big)
}
