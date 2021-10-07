use std::cmp::Ordering;

use crate::aabb::{surrounding_box, AABB};
use crate::hittable::{HitRecord, Hittable, HittableList};
use crate::ray::Ray;
use rand::Rng;


#[derive(Clone)]
pub struct BVHNode {
    pub left: Option<Box<dyn Hittable>>,
    pub right: Option<Box<dyn Hittable>>,
    pub aabb_box: AABB,
}

impl BVHNode {
    fn new(
        mut objects: Vec<Box<dyn Hittable>>,
        start: usize,
        end: usize,
        time0: f64,
        time1: f64,
    ) -> Self {
        let mut rnd = rand::thread_rng();

        let axis = rnd.gen_range(0, 3);
        let comparator = |x: &Box<dyn Hittable>, y: &Box<dyn Hittable>| {
            f64::partial_cmp(
                &(x.bounding_box(time0, time1).unwrap().min[axis]),
                &(y.bounding_box(time0, time1).unwrap().min[axis]),
            )
            .unwrap()
        };
        let left;
        let right;

        let object_span = end - start;

        if object_span == 1 {
            left = Some(objects[start].clone());
            right = Some(objects[start].clone());
        } else if object_span == 2 {
            let obj0 = objects[start].clone();
            let obj1 = objects[start + 1].clone();
            match comparator(&obj0, &obj1) {
                Ordering::Less => {
                    left = Some(obj0);
                    right = Some(obj1);
                }
                _ => {
                    left = Some(obj1);
                    right = Some(obj0);
                }
            }
        } else {
            objects.sort_unstable_by(comparator);
            let mid = start + object_span / 2;
            left = Some(Box::new(Self::new(objects.clone(), start, mid, time0, time1)));
            right = Some(Box::new(Self::new(objects.clone(), mid, end, time0, time1)));
        }

        if left.as_ref().unwrap().bounding_box(time0, time1).is_none()
            || right.as_ref().unwrap().bounding_box(time0, time1).is_none()
        {
            panic!("No bounding box in bvh_node constuctor.\n")
        }

        let aabb_box: AABB = surrounding_box(
            left.as_ref().unwrap().bounding_box(time0, time1).unwrap(),
            right.as_ref().unwrap().bounding_box(time0, time1).unwrap(),
        );

        return Self {
            left,
            right,
            aabb_box,
        };
    }
}

impl Hittable for BVHNode {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        if !self.aabb_box.hit(ray, t_min, t_max) {
            return None;
        }

        let mut hit_rec = None;
        let mut closest_so_far = t_max;

        if let Some(hit_left) = self.left.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
            closest_so_far = hit_left.t;
            hit_rec = Some(hit_left);
        }

        if self.right.is_some() {
            if let Some(hit_right) = self.right.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
                hit_rec = Some(hit_right);
            }
        }

        hit_rec
    }

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AABB> {
        Some(self.aabb_box)
    }
}
