use crate::aabb::AxisAlignedBoundingBox;
use crate::aarect::{XYRect, XZRect, YZRect};
use crate::hittable::{HitRecord, Hittable};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[allow(clippy::type_complexity)]
pub struct BoxEntity<M: Material + Clone> {
    box_min: Vec3,
    box_max: Vec3,
    sides: (
        XYRect<M>,
        XYRect<M>,
        XZRect<M>,
        XZRect<M>,
        YZRect<M>,
        YZRect<M>,
    ),
}

impl<M: Material + Clone> BoxEntity<M> {
    pub fn new(p0: Vec3, p1: Vec3, material: M) -> Self {
        Self {
            box_min: p0,
            box_max: p1,
            sides: (
                XYRect::new(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), material.clone()),
                XYRect::new(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), material.clone()),
                XZRect::new(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), material.clone()),
                XZRect::new(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), material.clone()),
                YZRect::new(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), material.clone()),
                YZRect::new(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), material),
            ),
        }
    }
}

macro_rules! hit {
    ($hitable:ident, $ray:ident, $t_min:ident, $closest_so_far:ident, $hit_record:ident) => {
        match $hitable.hit(&$ray, $t_min, $closest_so_far) {
            Some(hit_rec) => {
                $closest_so_far = hit_rec.t;
                $hit_record = Some(hit_rec);
            }
            None => {}
        }
    };
}

impl<M: Material + Clone> Hittable for BoxEntity<M> {
    #[allow(unused_assignments)]
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut hit_rec: Option<HitRecord> = None;
        let mut closest_so_far: f32 = t_max;

        let hitable = &self.sides.0;
        hit!(hitable, ray, t_min, closest_so_far, hit_rec);
        let hitable = &self.sides.1;
        hit!(hitable, ray, t_min, closest_so_far, hit_rec);
        let hitable = &self.sides.2;
        hit!(hitable, ray, t_min, closest_so_far, hit_rec);
        let hitable = &self.sides.3;
        hit!(hitable, ray, t_min, closest_so_far, hit_rec);
        let hitable = &self.sides.4;
        hit!(hitable, ray, t_min, closest_so_far, hit_rec);
        let hitable = &self.sides.5;
        hit!(hitable, ray, t_min, closest_so_far, hit_rec);
        hit_rec
    }

    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        Some(AxisAlignedBoundingBox::new(self.box_min, self.box_max))
    }
}
