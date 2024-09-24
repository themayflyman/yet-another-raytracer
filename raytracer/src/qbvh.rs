use std::simd::cmp::SimdPartialOrd;
use std::simd::num::SimdFloat;
use std::simd::{f32x4, u32x4};
use std::sync::Arc;

use crate::aabb::AxisAlignedBoundingBox;
use crate::hittable::{HitRecord, Hittable};
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct QBVH {
    hittables: Vec<Arc<dyn Hittable>>,
    nodes: Vec<QBVHNode>,
}

impl QBVH {
    pub fn new(hittables: Vec<Arc<dyn Hittable>>) -> Self {
        fn construct(
            hittables: &mut [Arc<dyn Hittable>],
            nodes: &mut Vec<QBVHNode>,
            indices: &mut [usize],
        ) -> (Option<AxisAlignedBoundingBox>, u32) {
            if hittables.len() == 0 {
                (None, u32::MAX)
            } else if hittables.len() == 1 {
                (
                    hittables[0].bounding_box(0.0, 0.0),
                    indices[0] as u32 | (1 << 31),
                )
            } else {
                let (left_hittables, left_indices, right_hittables, right_indices) =
                    split(hittables, indices);

                let (
                    left_left_hittables,
                    left_left_indices,
                    left_right_hittables,
                    left_right_indices,
                ) = split(left_hittables, left_indices);
                let (left_left_bbox, left_left_index) =
                    construct(left_left_hittables, nodes, left_left_indices);
                let (left_right_bbox, left_right_index) =
                    construct(left_right_hittables, nodes, left_right_indices);

                let (
                    right_left_hittables,
                    right_left_indices,
                    right_right_hittables,
                    right_right_indices,
                ) = split(right_hittables, right_indices);
                let (right_left_bbox, right_left_index) =
                    construct(right_left_hittables, nodes, right_left_indices);
                let (right_right_bbox, right_right_index) =
                    construct(right_right_hittables, nodes, right_right_indices);

                nodes.push(QBVHNode::new(
                    left_left_bbox,
                    left_left_index,
                    left_right_bbox,
                    left_right_index,
                    right_left_bbox,
                    right_left_index,
                    right_right_bbox,
                    right_right_index,
                ));

                let left_bbox = match (left_left_bbox, left_right_bbox) {
                    (Some(left_left_bbox), Some(left_right_bbox)) => {
                        AxisAlignedBoundingBox::surrounding_box(&left_left_bbox, &left_right_bbox)
                    }
                    (Some(left_left_bbox), None) => left_left_bbox,
                    (None, Some(left_right_bbox)) => left_right_bbox,
                    (None, None) => unreachable!(),
                };
                let right_bbox = match (right_left_bbox, right_right_bbox) {
                    (Some(right_left_bbox), Some(right_right_bbox)) => {
                        AxisAlignedBoundingBox::surrounding_box(&right_left_bbox, &right_right_bbox)
                    }
                    (Some(right_left_bbox), None) => right_left_bbox,
                    (None, Some(right_right_bbox)) => right_right_bbox,
                    (None, None) => unreachable!(),
                };

                (
                    Some(AxisAlignedBoundingBox::surrounding_box(
                        &left_bbox,
                        &right_bbox,
                    )),
                    (nodes.len() - 1) as u32,
                )
            }
        }

        let mut hittables: Vec<Arc<dyn Hittable>> = hittables.clone();
        let mut nodes: Vec<QBVHNode> = vec![];
        let mut indices: Vec<usize> = (0..hittables.len()).collect();

        construct(&mut hittables, &mut nodes, &mut indices);

        Self { hittables, nodes }
    }
}

impl Hittable for QBVH {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        Some(AxisAlignedBoundingBox::new(
            Vec3::new(
                self.nodes[0].bounding_box_min[0].reduce_min(),
                self.nodes[0].bounding_box_min[1].reduce_min(),
                self.nodes[0].bounding_box_min[2].reduce_min(),
            ),
            Vec3::new(
                self.nodes[0].bounding_box_max[0].reduce_max(),
                self.nodes[0].bounding_box_max[1].reduce_max(),
                self.nodes[0].bounding_box_max[2].reduce_max(),
            ),
        ))
    }

    fn hit(&self, ray: &crate::ray::Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        const MAX_STACK_SIZE: usize = 64;
        let nodes_len = self.nodes.len();
        let mut stack: [u32; MAX_STACK_SIZE] = [nodes_len as u32 - 1; MAX_STACK_SIZE];
        let mut stack_cursor: usize = 0;
        let mut hit_record: Option<HitRecord> = None;

        let ray_origin: [f32x4; 3] = [
            f32x4::splat(ray.origin().x()),
            f32x4::splat(ray.origin().y()),
            f32x4::splat(ray.origin().z()),
        ];
        let ray_direction: [f32x4; 3] = [
            f32x4::splat(ray.direction().x()),
            f32x4::splat(ray.direction().y()),
            f32x4::splat(ray.direction().z()),
        ];
        let inv_d: [f32x4; 3] = [
            f32x4::splat(1.0) / ray_direction[0],
            f32x4::splat(1.0) / ray_direction[1],
            f32x4::splat(1.0) / ray_direction[2],
        ];
        let t_minx4 = f32x4::splat(t_min);
        let mut t_max = t_max;

        loop {
            let id = stack[stack_cursor];
            if id >> 31 == 1 {
                let index = (id & ((1 << 31) - 1)) as usize;
                let hittiable = &self.hittables[index];
                if let Some(hr) = hittiable.hit(ray, t_min, t_max) {
                    t_max = t_max.min(hr.t);
                    hit_record = Some(hr);
                }
            } else {
                let node = &self.nodes[id as usize];
                let t_maxx4 = f32x4::splat(t_max);

                let t0: [f32x4; 3] = [
                    (node.bounding_box_min[0] - ray_origin[0]) * inv_d[0],
                    (node.bounding_box_min[1] - ray_origin[1]) * inv_d[1],
                    (node.bounding_box_min[2] - ray_origin[2]) * inv_d[2],
                ];
                let t1: [f32x4; 3] = [
                    (node.bounding_box_max[0] - ray_origin[0]) * inv_d[0],
                    (node.bounding_box_max[1] - ray_origin[1]) * inv_d[1],
                    (node.bounding_box_max[2] - ray_origin[2]) * inv_d[2],
                ];

                let t_minx4: f32x4 = [
                    t0[0].simd_min(t1[0]),
                    t0[1].simd_min(t1[1]),
                    t0[2].simd_min(t1[2]),
                ]
                .iter()
                .fold(t_minx4, |t_minx4, &t_curr| t_minx4.simd_max(t_curr));
                let t_maxx4: f32x4 = [
                    t0[0].simd_max(t1[0]),
                    t0[1].simd_max(t1[1]),
                    t0[2].simd_max(t1[2]),
                ]
                .iter()
                .fold(t_maxx4, |t_maxx4, &t_curr| t_maxx4.simd_min(t_curr));

                let hits = t_maxx4.simd_gt(t_minx4);
                for (index, hit) in (0..4).zip(hits.to_array().iter()) {
                    if *hit {
                        stack[stack_cursor] = node.children[index];
                        stack_cursor += 1;
                    }
                }
            }

            if stack_cursor == 0 {
                break;
            }
            stack_cursor -= 1;
        }

        hit_record
    }
}

#[derive(Clone)]
struct QBVHNode {
    bounding_box_min: [f32x4; 3],
    bounding_box_max: [f32x4; 3],
    children: u32x4, // sign is used to encode whether or not the child is a leaf node, the remaining 31 bits is used to encode the index to refer in triangles
}

impl QBVHNode {
    pub fn new(
        ll_bbox: Option<AxisAlignedBoundingBox>,
        ll_index: u32,
        lr_bbox: Option<AxisAlignedBoundingBox>,
        lr_index: u32,
        rl_bbox: Option<AxisAlignedBoundingBox>,
        rl_index: u32,
        rr_bbox: Option<AxisAlignedBoundingBox>,
        rr_index: u32,
    ) -> Self {
        let mut bounding_box_min = [f32x4::splat(f32::MAX); 3];
        let mut bounding_box_max = [f32x4::splat(f32::MAX); 3];
        let mut children = u32x4::splat(u32::MAX);

        let mut prepare_child_bbox_and_id =
            |bbox: Option<AxisAlignedBoundingBox>, offset: usize, index: u32| {
                if let Some(bbox) = bbox {
                    bounding_box_min[0][offset] = bbox.min[0];
                    bounding_box_min[1][offset] = bbox.min[1];
                    bounding_box_min[2][offset] = bbox.min[2];
                    bounding_box_max[0][offset] = bbox.max[0];
                    bounding_box_max[1][offset] = bbox.max[1];
                    bounding_box_max[2][offset] = bbox.max[2];
                }
                children[offset] = index;
            };
        prepare_child_bbox_and_id(ll_bbox, 0, ll_index);
        prepare_child_bbox_and_id(lr_bbox, 1, lr_index);
        prepare_child_bbox_and_id(rl_bbox, 2, rl_index);
        prepare_child_bbox_and_id(rr_bbox, 3, rr_index);

        Self {
            bounding_box_min,
            bounding_box_max,
            children,
        }
    }
}

#[inline(always)]
fn split<'a>(
    hittables: &'a mut [Arc<dyn Hittable>],
    indices: &'a mut [usize],
) -> (
    &'a mut [Arc<dyn Hittable>],
    &'a mut [usize],
    &'a mut [Arc<dyn Hittable>],
    &'a mut [usize],
) {
    let (min_x, max_x, min_y, max_y, min_z, max_z) = hittables.iter().fold(
        (
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::INFINITY,
            f32::NEG_INFINITY,
        ),
        |(min_x, max_x, min_y, max_y, min_z, max_z), hittable| {
            let centroid = hittable.centroid(0.0, 0.0).unwrap();
            (
                min_x.min(centroid.x()),
                max_x.max(centroid.x()),
                min_y.min(centroid.y()),
                max_y.max(centroid.y()),
                min_z.min(centroid.z()),
                max_z.max(centroid.z()),
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
    hittables.sort_unstable_by(|a, b| {
        f32::partial_cmp(
            &(a.centroid(0.0, 0.0).unwrap()[axis]),
            &(b.centroid(0.0, 0.0).unwrap()[axis]),
        )
        .unwrap()
    });

    let (hittables_left, hittables_right) = hittables.split_at_mut(hittables.len() / 2);
    let (indices_left, indices_right) = indices.split_at_mut(indices.len() / 2);

    (hittables_left, indices_left, hittables_right, indices_right)
}
