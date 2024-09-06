// use std::cmp::Ordering;

use crate::aabb::{surrounding_box, AxisAlignedBoundingBox};
use crate::hittable::{HitRecord, Hittable};
use crate::ray::Ray;
// use rand::Rng;
use std::sync::Arc;

#[derive(Clone)]
pub struct BVHNode {
    pub left: Option<Arc<dyn Hittable>>,
    pub right: Option<Arc<dyn Hittable>>,
    pub aabb_box: AxisAlignedBoundingBox,
}

impl BVHNode {
    pub fn new(
        objects: &mut Vec<Arc<dyn Hittable>>,
        start: usize,
        end: usize,
        time0: f32,
        time1: f32,
    ) -> Self {
        // let mut rnd = rand::thread_rng();

        // let axis = rnd.gen_range(0, 3);
        // let comparator = |a: &Arc<dyn Hittable>, b: &Arc<dyn Hittable>| {
        //     let max_x = [a, b].iter().fold(f32::NEG_INFINITY, |max, &v| {
        //         let centroid = v.centroid(time0, time1).unwrap();
        //         return f32::max(max, centroid.x());
        //     });
        //     let max_y = [a, b].iter().fold(f32::NEG_INFINITY, |max, &v| {
        //         let centroid = v.centroid(time0, time1).unwrap();
        //         return f32::max(max, centroid.y());
        //     });
        //     let max_z = [a, b].iter().fold(f32::NEG_INFINITY, |max, &v| {
        //         let centroid = v.centroid(time0, time1).unwrap();
        //         return f32::max(max, centroid.z());
        //     });
        //     let min_x = [a, b].iter().fold(f32::INFINITY, |min, &v| {
        //         let centroid = v.centroid(time0, time1).unwrap();
        //         return f32::min(min, centroid.x());
        //     });
        //     let min_y = [a, b].iter().fold(f32::INFINITY, |min, &v| {
        //         let centroid = v.centroid(time0, time1).unwrap();
        //         return f32::min(min, centroid.y());
        //     });
        //     let min_z = [a, b].iter().fold(f32::INFINITY, |min, &v| {
        //         let centroid = v.centroid(time0, time1).unwrap();
        //         return f32::min(min, centroid.z());
        //     });

        //     let mut axis: usize = 0;
        //     if max_y - min_y > max_x - min_x {
        //         axis = 1;
        //     }
        //     if max_z - min_z > f32::max(max_y - min_y, max_x - min_x) {
        //         axis = 2;
        //     }

        //     f32::partial_cmp(
        //         &(a.centroid(time0, time1).unwrap()[axis]),
        //         &(b.centroid(time0, time1).unwrap()[axis]),
        //     )
        //     .unwrap()
        // };
        let left;
        let right;

        let (min_x, max_x, min_y, max_y, min_z, max_z) = objects.iter().fold(
            (
                f32::INFINITY,
                f32::NEG_INFINITY,
                f32::INFINITY,
                f32::NEG_INFINITY,
                f32::INFINITY,
                f32::NEG_INFINITY,
            ),
            |(min_x, max_x, min_y, max_y, min_z, max_z), obj| {
                let centroid = obj.centroid(time0, time1).unwrap();
                (
                    f32::min(min_x, centroid.x()),
                    f32::max(max_x, centroid.x()),
                    f32::min(min_y, centroid.y()),
                    f32::max(max_y, centroid.y()),
                    f32::min(min_z, centroid.z()),
                    f32::max(max_z, centroid.z()),
                )
            },
        );
        let mut axis: usize = 0;
        if max_y - min_y > max_x - min_x {
            axis = 1;
        }
        if max_z - min_z > f32::max(max_y - min_y, max_x - min_x) {
            axis = 2;
        }
        objects.sort_unstable_by(|a, b| {
            f32::partial_cmp(
                &(a.centroid(time0, time1).unwrap()[axis]),
                &(b.centroid(time0, time1).unwrap()[axis]),
            )
            .unwrap()
        });

        let object_span = end - start;
        if object_span == 1 {
            left = Some(objects[start].clone());
            right = Some(objects[start].clone());
        } else if object_span == 2 {
            // let obj0 = objects[start].clone();
            // let obj1 = objects[start + 1].clone();
            // match comparator(&obj0, &obj1) {
            //     Ordering::Less => {
            //         left = Some(obj0);
            //         right = Some(obj1);
            //     }
            //     _ => {
            //         left = Some(obj1);
            //         right = Some(obj0);
            //     }
            // }
            left = Some(objects[start].clone());
            right = Some(objects[start + 1].clone());
        } else {
            // objects.sort_unstable_by(comparator);
            let mid = start + object_span / 2;
            left = Some(Arc::new(Self::new(objects, start, mid, time0, time1)));
            right = Some(Arc::new(Self::new(objects, mid, end, time0, time1)));
        }

        if left.as_ref().unwrap().bounding_box(time0, time1).is_none()
            || right.as_ref().unwrap().bounding_box(time0, time1).is_none()
        {
            panic!("No bounding box in bvh_node constuctor.\n")
        }

        let aabb_box: AxisAlignedBoundingBox = surrounding_box(
            left.as_ref().unwrap().bounding_box(time0, time1).unwrap(),
            right.as_ref().unwrap().bounding_box(time0, time1).unwrap(),
        );

        Self {
            left,
            right,
            aabb_box,
        }
    }
}

impl Hittable for BVHNode {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        // let mut hit_rec = None;
        // let mut closest_so_far = t_max;

        // if let Some(hit_left) = self.left.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
        //     closest_so_far = hit_left.t;
        //     hit_rec = Some(hit_left);
        // }

        // if let Some(hit_right) = self.right.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
        //     hit_rec = Some(hit_right);
        // }

        // hit_rec

        let left = self.left.as_ref().unwrap();
        let right = self.right.as_ref().unwrap();

        let (hit_left, hit_right) = left.bounding_box(0.0, 1.0).unwrap().hit_with_another(
            &right.bounding_box(0.0, 1.0).unwrap(),
            ray,
            t_min,
            t_max,
        );

        match(hit_left, hit_right) {
            (None, None) => None,
            (Some(_t), None) => left.hit(ray, t_min, t_max),
            (None, Some(_t)) => right.hit(ray, t_min, t_max),
            (Some(t_left), Some(t_right)) => {
                let mut hit_rec = None;
                let mut closest_so_far = t_max;

                if t_left < t_right {
                  if let Some(hit_left) = self.left.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
                      closest_so_far = hit_left.t;
                      hit_rec = Some(hit_left);
                  }

                  if let Some(hit_right) = self.right.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
                      hit_rec = Some(hit_right);
                  }
                } else {
                  if let Some(hit_right) = self.right.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
                      closest_so_far = hit_right.t;
                      hit_rec = Some(hit_right);
                  }

                  if let Some(hit_left) = self.left.as_ref().unwrap().hit(ray, t_min, closest_so_far) {
                      hit_rec = Some(hit_left);
                  }
                }
                hit_rec
            },
        }
    }

    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        Some(self.aabb_box)
    }
}

pub struct BVHNodeStatic<L: Hittable, R: Hittable> {
    pub left: Arc<L>,
    pub right: Arc<R>,
    pub aabb_box: AxisAlignedBoundingBox,
}

impl<L: Hittable, R: Hittable> BVHNodeStatic<L, R> {
    pub fn construct(left: Arc<L>, right: Arc<R>) -> Self {
        let aabb_box = AxisAlignedBoundingBox::surrounding_box(
            left.bounding_box(0.0, 1.0).as_ref().unwrap(),
            right.bounding_box(0.0, 1.0).as_ref().unwrap(),
        );
        Self {
            left,
            right,
            aabb_box,
        }
    }
}

impl<L: Hittable, R: Hittable> Hittable for BVHNodeStatic<L, R> {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        match self.aabb_box.hit(ray, t_min, t_max) {
            false => None,
            true => {
                let hit_left = self.left.hit(ray, t_min, t_max);
                let hit_right = self.right.hit(ray, t_min, t_max);
                match (hit_left, hit_right) {
                    (None, None) => None,
                    (None, Some(hit_rec)) => Some(hit_rec),
                    (Some(hit_rec), None) => Some(hit_rec),
                    (Some(hit_left), Some(hit_right)) => {
                        if hit_left.t < hit_right.t {
                            Some(hit_left)
                        } else {
                            Some(hit_right)
                        }
                    }
                }
            }
        }
    }

    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        Some(self.aabb_box)
    }
}
