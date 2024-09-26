use std::collections::HashMap;
use std::ops::Not;
use std::simd::cmp::SimdPartialOrd;
use std::simd::num::SimdFloat;
use std::simd::{f32x4, u32x4};
use std::sync::Arc;

use crate::aabb::AxisAlignedBoundingBox;
use crate::hittable::{HitRecord, Hittable};
use crate::material::Material;
use crate::triangle::Triangle;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct L1QBVH<M: Material> {
    triangles: Vec<Arc<Triangle<M>>>,
    nodes: Vec<QBVHNode>,
}

impl<M: Material> L1QBVH<M> {
    pub fn new(triangles: Vec<Arc<Triangle<M>>>) -> Self {
        fn construct<M: Material>(
            triangles: &mut [Arc<Triangle<M>>],
            nodes: &mut Vec<QBVHNode>,
            indices: &mut [usize],
        ) -> (Option<AxisAlignedBoundingBox>, u32) {
            if triangles.len() == 0 {
                (None, u32::MAX)
            } else if triangles.len() == 1 {
                (
                    triangles[0].bounding_box(0.0, 0.0),
                    indices[0] as u32 | (1 << 31),
                )
            } else {
                let (left_triangles, left_indices, right_triangles, right_indices) =
                    split(triangles, indices);

                let (
                    left_left_triangles,
                    left_left_indices,
                    left_right_triangles,
                    left_right_indices,
                ) = split(left_triangles, left_indices);
                let (left_left_bbox, left_left_index) =
                    construct(left_left_triangles, nodes, left_left_indices);
                let (left_right_bbox, left_right_index) =
                    construct(left_right_triangles, nodes, left_right_indices);

                let (
                    right_left_triangles,
                    right_left_indices,
                    right_right_triangles,
                    right_right_indices,
                ) = split(right_triangles, right_indices);
                let (right_left_bbox, right_left_index) =
                    construct(right_left_triangles, nodes, right_left_indices);
                let (right_right_bbox, right_right_index) =
                    construct(right_right_triangles, nodes, right_right_indices);

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

        let mut triangles: Vec<Arc<Triangle<M>>> = triangles.clone();
        let mut nodes: Vec<QBVHNode> = vec![];
        let mut indices: Vec<usize> = (0..triangles.len()).collect();

        construct(&mut triangles, &mut nodes, &mut indices);

        Self {
            triangles,
            nodes,
        }
    }
}

impl<M: Material> Hittable for L1QBVH<M> {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        let nodes_len = self.nodes.len();
        Some(AxisAlignedBoundingBox::new(
            Vec3::new(
                self.nodes[nodes_len - 1].bounding_box_min[0].reduce_min(),
                self.nodes[nodes_len - 1].bounding_box_min[1].reduce_min(),
                self.nodes[nodes_len - 1].bounding_box_min[2].reduce_min(),
            ),
            Vec3::new(
                self.nodes[nodes_len - 1].bounding_box_max[0].reduce_max(),
                self.nodes[nodes_len - 1].bounding_box_max[1].reduce_max(),
                self.nodes[nodes_len - 1].bounding_box_max[2].reduce_max(),
            ),
        ))
    }

    fn hit(&self, ray: &crate::ray::Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        const MAX_STACK_SIZE: usize = 64;
        let nodes_len = self.nodes.len();
        let mut stack: [u32; MAX_STACK_SIZE] = [nodes_len as u32 - 1; MAX_STACK_SIZE];
        let mut stack_cursor: usize = 0;
        let mut hit_record: Option<HitRecord> = None;

        let rox4: [f32x4; 3] = [
            f32x4::splat(ray.origin().x()),
            f32x4::splat(ray.origin().y()),
            f32x4::splat(ray.origin().z()),
        ];
        let rdx4: [f32x4; 3] = [
            f32x4::splat(ray.direction().x()),
            f32x4::splat(ray.direction().y()),
            f32x4::splat(ray.direction().z()),
        ];
        let inv_rdx4: [f32x4; 3] = [
            f32x4::splat(1.0) / rdx4[0],
            f32x4::splat(1.0) / rdx4[1],
            f32x4::splat(1.0) / rdx4[2],
        ];
        let t_minx4 = f32x4::splat(t_min);
        let mut t_max = t_max;

        loop {
            let id = stack[stack_cursor];
            if id >> 31 == 1 {
                let index = (id & ((1 << 31) - 1)) as usize;
                let hittiable = &self.triangles[index];
                if let Some(hr) = hittiable.hit(ray, t_min, t_max) {
                    t_max = t_max.min(hr.t);
                    hit_record = Some(hr);
                }
            } else {
                let node = &self.nodes[id as usize];
                let t_maxx4 = f32x4::splat(t_max);

                let t0: [f32x4; 3] = [
                    (node.bounding_box_min[0] - rox4[0]) * inv_rdx4[0],
                    (node.bounding_box_min[1] - rox4[1]) * inv_rdx4[1],
                    (node.bounding_box_min[2] - rox4[2]) * inv_rdx4[2],
                ];
                let t1: [f32x4; 3] = [
                    (node.bounding_box_max[0] - rox4[0]) * inv_rdx4[0],
                    (node.bounding_box_max[1] - rox4[1]) * inv_rdx4[1],
                    (node.bounding_box_max[2] - rox4[2]) * inv_rdx4[2],
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

                let hits = t_maxx4.simd_gt(t_minx4).to_int();
                for i in 0..4 {
                    stack[stack_cursor] = node.children[i];
                    stack_cursor -= hits[i] as usize;
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
pub struct L4QBVH<M: Material> {
    triangles: Vec<Arc<Triangle<M>>>,
    nodes: Vec<QBVHNode>,
    soa: HashMap<u32, [[f32x4; 3]; 8]>,
}

impl<M: Material> L4QBVH<M> {
    pub fn new(triangles: Vec<Arc<Triangle<M>>>) -> Self {
        fn construct<M: Material>(
            triangles: &mut [Arc<Triangle<M>>],
            nodes: &mut Vec<QBVHNode>,
            indices: &mut [usize],
            soa: &mut HashMap<u32, [[f32x4; 3]; 8]>,
        ) -> (Option<AxisAlignedBoundingBox>, u32) {
            if triangles.len() == 0 {
                (None, u32::MAX)
            } else if triangles.len() <= 4 {
                let bounding_box = triangles.iter().fold(None, |bbox, triangle| match bbox {
                    Some(bbox) => Some(AxisAlignedBoundingBox::surrounding_box(
                        &bbox,
                        &triangle.bounding_box(0.0, 0.0).unwrap(),
                    )),
                    None => Some(triangle.bounding_box(0.0, 0.0).unwrap()),
                });

                let id = indices[0] as u32 | (1 << 31) | ((triangles.len() as u32) << 27);
                let mut t_v0x4 = [f32x4::splat(f32::MAX); 3];
                let mut t_v1x4 = [f32x4::splat(f32::MAX); 3];
                let mut t_v2x4 = [f32x4::splat(f32::MAX); 3];
                let mut t_normals0x4 = [f32x4::splat(f32::MAX); 3];
                let mut t_normals1x4 = [f32x4::splat(f32::MAX); 3];
                let mut t_normals2x4 = [f32x4::splat(f32::MAX); 3];
                let mut t_ux4 = [f32x4::splat(f32::MAX); 3];
                let mut t_vx4 = [f32x4::splat(f32::MAX); 3];
                for (i, triangle) in triangles.iter().enumerate() {
                    for j in 0..3 {
                        t_v0x4[j][i] = triangle.vertices[0][j];
                        t_v1x4[j][i] = triangle.vertices[1][j];
                        t_v2x4[j][i] = triangle.vertices[2][j];
                        t_normals0x4[j][i] = triangle.normals[0][j];
                        t_normals1x4[j][i] = triangle.normals[1][j];
                        t_normals2x4[j][i] = triangle.normals[2][j];
                        t_ux4[j][i] = triangle.uv[j].0;
                        t_vx4[j][i] = triangle.uv[j].1;
                    }
                }
                let edge1x4 = [
                    t_v1x4[0] - t_v0x4[0],
                    t_v1x4[1] - t_v0x4[1],
                    t_v1x4[2] - t_v0x4[2],
                ];
                let edge2x4 = [
                    t_v2x4[0] - t_v0x4[0],
                    t_v2x4[1] - t_v0x4[1],
                    t_v2x4[2] - t_v0x4[2],
                ];
                soa.insert(
                    id,
                    [
                        t_v0x4,
                        edge1x4,
                        edge2x4,
                        t_normals0x4,
                        t_normals1x4,
                        t_normals2x4,
                        t_ux4,
                        t_vx4,
                    ],
                );

                (
                    bounding_box,
                    // encode high bit to indicate that this is a leaf node, use the next 4 bit
                    // encode the len of the triangles, use the rest of the bits to encode the
                    // indices[0]
                    id,
                )
            } else {
                let (left_triangles, left_indices, right_triangles, right_indices) =
                    split(triangles, indices);

                let (
                    left_left_triangles,
                    left_left_indices,
                    left_right_triangles,
                    left_right_indices,
                ) = split(left_triangles, left_indices);
                let (left_left_bbox, left_left_index) =
                    construct(left_left_triangles, nodes, left_left_indices, soa);
                let (left_right_bbox, left_right_index) =
                    construct(left_right_triangles, nodes, left_right_indices, soa);

                let (
                    right_left_triangles,
                    right_left_indices,
                    right_right_triangles,
                    right_right_indices,
                ) = split(right_triangles, right_indices);
                let (right_left_bbox, right_left_index) =
                    construct(right_left_triangles, nodes, right_left_indices, soa);
                let (right_right_bbox, right_right_index) =
                    construct(right_right_triangles, nodes, right_right_indices, soa);

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

        let mut triangles: Vec<Arc<Triangle<M>>> = triangles.clone();
        let mut nodes: Vec<QBVHNode> = vec![];
        let mut indices: Vec<usize> = (0..triangles.len()).collect();
        let mut soa: HashMap<u32, [[f32x4; 3]; 8]> = HashMap::new();

        construct(&mut triangles, &mut nodes, &mut indices, &mut soa);

        Self {
            triangles,
            nodes,
            soa,
        }
    }
}

impl<M: Material> Hittable for L4QBVH<M> {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        let nodes_len = self.nodes.len();
        Some(AxisAlignedBoundingBox::new(
            Vec3::new(
                self.nodes[nodes_len - 1].bounding_box_min[0].reduce_min(),
                self.nodes[nodes_len - 1].bounding_box_min[1].reduce_min(),
                self.nodes[nodes_len - 1].bounding_box_min[2].reduce_min(),
            ),
            Vec3::new(
                self.nodes[nodes_len - 1].bounding_box_max[0].reduce_max(),
                self.nodes[nodes_len - 1].bounding_box_max[1].reduce_max(),
                self.nodes[nodes_len - 1].bounding_box_max[2].reduce_max(),
            ),
        ))
    }

    fn hit(&self, ray: &crate::ray::Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        const MAX_STACK_SIZE: usize = 64;
        let nodes_len = self.nodes.len();
        let mut stack: [u32; MAX_STACK_SIZE] = [nodes_len as u32 - 1; MAX_STACK_SIZE];
        let mut stack_cursor: usize = 0;
        let mut hit_record: Option<HitRecord> = None;

        let rox4: [f32x4; 3] = [
            f32x4::splat(ray.origin().x()),
            f32x4::splat(ray.origin().y()),
            f32x4::splat(ray.origin().z()),
        ];
        let rdx4: [f32x4; 3] = [
            f32x4::splat(ray.direction().x()),
            f32x4::splat(ray.direction().y()),
            f32x4::splat(ray.direction().z()),
        ];
        let inv_rdx4: [f32x4; 3] = [
            f32x4::splat(1.0) / rdx4[0],
            f32x4::splat(1.0) / rdx4[1],
            f32x4::splat(1.0) / rdx4[2],
        ];
        let t_minx4 = f32x4::splat(t_min);
        let mut t_max = t_max;

        loop {
            let id = stack[stack_cursor];
            if id >> 31 == 1 {
                let count = ((id & (0b1111 << 27)) >> 27) as usize;
                let index = (id & ((1 << 27) - 1)) as usize;

                let [t_v0x4, edge1x4, edge2x4, t_normals0x4, t_normals1x4, t_normals2x4, t_ux4, t_vx4] =
                    self.soa.get(&(id as u32)).unwrap();

                let hx4 = [
                    rdx4[1] * edge2x4[2] - rdx4[2] * edge2x4[1],
                    rdx4[2] * edge2x4[0] - rdx4[0] * edge2x4[2],
                    rdx4[0] * edge2x4[1] - rdx4[1] * edge2x4[0],
                ];
                let ax4 = edge1x4[0] * hx4[0] + edge1x4[1] * hx4[1] + edge1x4[2] * hx4[2];

                let mut hitx4 = (ax4.simd_gt(-f32x4::splat(f32::EPSILON))
                    & ax4.simd_lt(f32x4::splat(f32::EPSILON)))
                .not();

                let fx4 = f32x4::splat(1.0) / ax4;
                let sx4 = [
                    rox4[0] - t_v0x4[0],
                    rox4[1] - t_v0x4[1],
                    rox4[2] - t_v0x4[2],
                ];
                let ux4 = fx4 * (sx4[0] * hx4[0] + sx4[1] * hx4[1] + sx4[2] * hx4[2]);

                hitx4 &= ux4.simd_ge(f32x4::splat(0.0)) & ux4.simd_le(f32x4::splat(1.0));

                let qx4 = [
                    sx4[1] * edge1x4[2] - sx4[2] * edge1x4[1],
                    sx4[2] * edge1x4[0] - sx4[0] * edge1x4[2],
                    sx4[0] * edge1x4[1] - sx4[1] * edge1x4[0],
                ];
                let vx4 = fx4 * (rdx4[0] * qx4[0] + rdx4[1] * qx4[1] + rdx4[2] * qx4[2]);
                hitx4 &= vx4.simd_ge(f32x4::splat(0.0)) & (ux4 + vx4).simd_le(f32x4::splat(1.0));

                let tx4 = fx4 * (edge2x4[0] * qx4[0] + edge2x4[1] * qx4[1] + edge2x4[2] * qx4[2]);
                hitx4 &= tx4.simd_ge(t_minx4) & tx4.simd_le(f32x4::splat(t_max));

                let px4 = [
                    rox4[0] + tx4 * rdx4[0],
                    rox4[1] + tx4 * rdx4[1],
                    rox4[2] + tx4 * rdx4[2],
                ];
                let wx4 = f32x4::splat(1.0) - ux4 - vx4;
                let outward_normalx4 = [
                    t_normals0x4[0] * wx4 + t_normals1x4[0] * ux4 + t_normals2x4[0] * vx4,
                    t_normals0x4[1] * wx4 + t_normals1x4[1] * ux4 + t_normals2x4[1] * vx4,
                    t_normals0x4[2] * wx4 + t_normals1x4[2] * ux4 + t_normals2x4[2] * vx4,
                ];
                let ffx4 = (rdx4[0] * outward_normalx4[0]
                    + rdx4[1] * outward_normalx4[1]
                    + rdx4[2] * outward_normalx4[2])
                    .simd_le(f32x4::splat(0.0));
                let h_ux4 = t_ux4[0] * wx4 + t_ux4[1] * ux4 + t_ux4[2] * vx4;
                let h_vx4 = t_vx4[0] * wx4 + t_vx4[1] * ux4 + t_vx4[2] * vx4;

                for i in 0..count {
                    let ff = ffx4.test(i);
                    let sign: f32 = if ff { 1.0 } else { -1.0 };
                    let normal = Vec3::new(
                        sign * outward_normalx4[0][i],
                        sign * outward_normalx4[1][i],
                        sign * outward_normalx4[2][i],
                    );
                    if hitx4.test(i) && t_max > tx4[i] {
                        t_max = tx4[i];
                        hit_record = Some(HitRecord {
                            u: h_ux4[i],
                            v: h_vx4[i],
                            t: tx4[i],
                            p: Vec3::new(px4[0][i], px4[1][i], px4[2][i]),
                            normal,
                            front_face: ff,
                            material: &self.triangles[index + i].material,
                        })
                    }
                }
            } else {
                let node = &self.nodes[id as usize];
                let t_maxx4 = f32x4::splat(t_max);

                let t0: [f32x4; 3] = [
                    (node.bounding_box_min[0] - rox4[0]) * inv_rdx4[0],
                    (node.bounding_box_min[1] - rox4[1]) * inv_rdx4[1],
                    (node.bounding_box_min[2] - rox4[2]) * inv_rdx4[2],
                ];
                let t1: [f32x4; 3] = [
                    (node.bounding_box_max[0] - rox4[0]) * inv_rdx4[0],
                    (node.bounding_box_max[1] - rox4[1]) * inv_rdx4[1],
                    (node.bounding_box_max[2] - rox4[2]) * inv_rdx4[2],
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

                let hits = t_maxx4.simd_gt(t_minx4).to_int();
                for i in 0..4 {
                    stack[stack_cursor] = node.children[i];
                    stack_cursor -= hits[i] as usize;
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
fn split<'a, M: Material>(
    triangles: &'a mut [Arc<Triangle<M>>],
    indices: &'a mut [usize],
) -> (
    &'a mut [Arc<Triangle<M>>],
    &'a mut [usize],
    &'a mut [Arc<Triangle<M>>],
    &'a mut [usize],
) {
    let (min_x, max_x, min_y, max_y, min_z, max_z) = triangles.iter().fold(
        (
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::INFINITY,
            f32::NEG_INFINITY,
        ),
        |(min_x, max_x, min_y, max_y, min_z, max_z), triangle| {
            let centroid = triangle.centroid(0.0, 0.0).unwrap();
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
    triangles.sort_unstable_by(|a, b| {
        f32::partial_cmp(
            &(a.centroid(0.0, 0.0).unwrap()[axis]),
            &(b.centroid(0.0, 0.0).unwrap()[axis]),
        )
        .unwrap()
    });

    let (triangles_left, triangles_right) = triangles.split_at_mut(triangles.len() / 2);
    let (indices_left, indices_right) = indices.split_at_mut(indices.len() / 2);

    (triangles_left, indices_left, triangles_right, indices_right)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bvh::BVHNode;
    use crate::color::RGB;
    use crate::material::Lambertian;
    use crate::ray::Ray;
    use crate::texture::SolidColor;
    use crate::vec3::Vec3;
    use rand::prelude::*;
    use std::path::Path;
    use test::{black_box, Bencher};

    fn load_hittables_from_obj(obj_path: &Path, material: impl Clone + Material + 'static) -> Vec<Arc<dyn Hittable>> {
        let (models, _materials) =
            tobj::load_obj(obj_path, &tobj::GPU_LOAD_OPTIONS).expect("Failed to load OBJ file");

        let mut triangles: Vec<Arc<dyn Hittable>> = vec![];
        for model in models.iter() {
            let mesh = &model.mesh;
            let mesh_vertices: Vec<Vec3> = mesh
                .positions
                .chunks_exact(3)
                .map(|p| Vec3::new(p[0].into(), p[1].into(), p[2].into()))
                .collect();
            let mesh_normals: Vec<Option<Vec3>> = if !mesh.normals.is_empty() {
                mesh.normals
                    .chunks_exact(3)
                    .map(|n| Some(Vec3::new(n[0].into(), n[1].into(), n[2].into())))
                    .collect()
            } else {
                mesh.indices.iter().map(|_| None).collect()
            };
            let mesh_uv: Vec<Option<(f32, f32)>> = if !mesh.texcoords.is_empty() {
                mesh.texcoords
                    .chunks_exact(2)
                    .map(|uv| Some((uv[0].into(), uv[1].into())))
                    .collect()
            } else {
                mesh.indices.iter().map(|_| None).collect()
            };

            for i in 0..mesh.indices.len() / 3 {
                let vertices = [
                    mesh_vertices[mesh.indices[i * 3] as usize],
                    mesh_vertices[mesh.indices[i * 3 + 1] as usize],
                    mesh_vertices[mesh.indices[i * 3 + 2] as usize],
                ];
                let default_normal = (vertices[1] - vertices[0])
                    .cross(&(vertices[2] - vertices[0]))
                    .unit_vector();
                let normals = [
                    mesh_normals[mesh.indices[i * 3] as usize].unwrap_or(default_normal),
                    mesh_normals[mesh.indices[i * 3 + 1] as usize].unwrap_or(default_normal),
                    mesh_normals[mesh.indices[i * 3 + 2] as usize].unwrap_or(default_normal),
                ];
                let uv = [
                    mesh_uv[mesh.indices[i * 3] as usize].unwrap_or((0.0, 0.0)),
                    mesh_uv[mesh.indices[i * 3 + 1] as usize].unwrap_or((0.0, 0.0)),
                    mesh_uv[mesh.indices[i * 3 + 2] as usize].unwrap_or((0.0, 0.0)),
                ];

                let triangle = Triangle {
                    vertices,
                    normals,
                    uv,
                    material: material.clone(),
                };
                triangles.push(Arc::new(triangle));
            }
        }

        triangles
    }

    fn load_triangles_from_obj<M: Clone + Material>(obj_path: &Path, material: M) -> Vec<Arc<Triangle<M>>> {
        let (models, _materials) =
            tobj::load_obj(obj_path, &tobj::GPU_LOAD_OPTIONS).expect("Failed to load OBJ file");

        let mut triangles: Vec<Arc<Triangle<M>>> = vec![];
        for model in models.iter() {
            let mesh = &model.mesh;
            let mesh_vertices: Vec<Vec3> = mesh
                .positions
                .chunks_exact(3)
                .map(|p| Vec3::new(p[0].into(), p[1].into(), p[2].into()))
                .collect();
            let mesh_normals: Vec<Option<Vec3>> = if !mesh.normals.is_empty() {
                mesh.normals
                    .chunks_exact(3)
                    .map(|n| Some(Vec3::new(n[0].into(), n[1].into(), n[2].into())))
                    .collect()
            } else {
                mesh.indices.iter().map(|_| None).collect()
            };
            let mesh_uv: Vec<Option<(f32, f32)>> = if !mesh.texcoords.is_empty() {
                mesh.texcoords
                    .chunks_exact(2)
                    .map(|uv| Some((uv[0].into(), uv[1].into())))
                    .collect()
            } else {
                mesh.indices.iter().map(|_| None).collect()
            };

            for i in 0..mesh.indices.len() / 3 {
                let vertices = [
                    mesh_vertices[mesh.indices[i * 3] as usize],
                    mesh_vertices[mesh.indices[i * 3 + 1] as usize],
                    mesh_vertices[mesh.indices[i * 3 + 2] as usize],
                ];
                let default_normal = (vertices[1] - vertices[0])
                    .cross(&(vertices[2] - vertices[0]))
                    .unit_vector();
                let normals = [
                    mesh_normals[mesh.indices[i * 3] as usize].unwrap_or(default_normal),
                    mesh_normals[mesh.indices[i * 3 + 1] as usize].unwrap_or(default_normal),
                    mesh_normals[mesh.indices[i * 3 + 2] as usize].unwrap_or(default_normal),
                ];
                let uv = [
                    mesh_uv[mesh.indices[i * 3] as usize].unwrap_or((0.0, 0.0)),
                    mesh_uv[mesh.indices[i * 3 + 1] as usize].unwrap_or((0.0, 0.0)),
                    mesh_uv[mesh.indices[i * 3 + 2] as usize].unwrap_or((0.0, 0.0)),
                ];

                let triangle = Triangle {
                    vertices,
                    normals,
                    uv,
                    material: material.clone(),
                };
                triangles.push(Arc::new(triangle));
            }
        }

        triangles
    }

    fn generate_worse_case_rays(
        count: usize,
        bounds: &AxisAlignedBoundingBox,
    ) -> Vec<crate::ray::Ray> {
        let mut rng = rand::thread_rng();
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];

        (0..count)
            .map(|_| {
                let origin = Vec3::new(
                    rng.gen_range(bounds.min.x()..bounds.max.x()),
                    rng.gen_range(bounds.min.y()..bounds.max.y()),
                    rng.gen_range(bounds.min.z()..bounds.max.z()),
                );
                let direction = axes[rng.gen_range(0..3)];
                Ray::new(origin, direction, 0.0, 0.0)
            })
            .collect()
    }

    fn generate_uniform_random_rays(
        count: usize,
        bounds: &AxisAlignedBoundingBox
    ) -> Vec<Ray> {
        let mut rng = rand::thread_rng();
        (0..count)
            .map(|_| {
                let origin = Vec3::new(
                    rng.gen_range(bounds.min.x()..bounds.max.x()),
                    rng.gen_range(bounds.min.y()..bounds.max.y()),
                    rng.gen_range(bounds.min.z()..bounds.max.z()),
                );
                let direction = Vec3::random_unit_vector(&mut rng);
                Ray::new(origin, direction, 0.0, 0.0)
            })
            .collect()
    }

    #[bench]
    fn bench_hit_l1_qbvh_stanford_bunny_worse_case_rays(bencher: &mut Bencher) {
        let material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
        let stanford_bunny = load_triangles_from_obj(Path::new("../input/bunny.obj"), material);
        let l1_qbvh = L1QBVH::new(stanford_bunny);

        let rays = generate_worse_case_rays(512, &l1_qbvh.bounding_box(0.0, 0.0).unwrap());
        let mut rays_in = rays.iter().cycle();

        bencher.iter(|| {
            let ray = rays_in.next().unwrap();
            black_box(l1_qbvh.hit(&ray, 0.0, f32::INFINITY));
        });
    }

    #[bench]
    fn bench_hit_l4_qbvh_stanford_bunny_worse_case_rays(bencher: &mut Bencher) {
        let material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
        let stanford_bunny = load_triangles_from_obj(Path::new("../input/bunny.obj"), material);
        let l4_qbvh = L4QBVH::new(stanford_bunny);

        let rays = generate_worse_case_rays(512, &l4_qbvh.bounding_box(0.0, 0.0).unwrap());
        let mut rays_in = rays.iter().cycle();

        bencher.iter(|| {
            let ray = rays_in.next().unwrap();
            black_box(l4_qbvh.hit(&ray, 0.0, f32::INFINITY));
        });
    }

    #[bench]
    fn bench_hit_l1_qbvh_stanford_bunny_uniform_rays(bencher: &mut Bencher) {
        let material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
        let stanford_bunny = load_triangles_from_obj(Path::new("../input/bunny.obj"), material);
        let l1_qbvh = L1QBVH::new(stanford_bunny);

        let rays = generate_uniform_random_rays(512, &l1_qbvh.bounding_box(0.0, 0.0).unwrap());
        let mut rays_in = rays.iter().cycle();

        bencher.iter(|| {
            let ray = rays_in.next().unwrap();
            black_box(l1_qbvh.hit(&ray, 0.0, f32::INFINITY));
        });
    }

    #[bench]
    fn bench_hit_l4_qbvh_stanford_bunny_uniform_rays(bencher: &mut Bencher) {
        let material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
        let stanford_bunny = load_triangles_from_obj(Path::new("../input/bunny.obj"), material);
        let l4_qbvh = L4QBVH::new(stanford_bunny);

        let rays = generate_uniform_random_rays(512, &l4_qbvh.bounding_box(0.0, 0.0).unwrap());
        let mut rays_in = rays.iter().cycle();

        bencher.iter(|| {
            let ray = rays_in.next().unwrap();
            black_box(l4_qbvh.hit(&ray, 0.0, f32::INFINITY));
        });
    }

    // #[bench]
    // fn bench_hit_bvh_stanford_bunny_worse_case_rays(bencher: &mut Bencher) {
    //     let material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
    //     let mut stanford_bunny = load_hittables_from_obj(Path::new("../input/bunny.obj"), material);
    //     let size = stanford_bunny.len();
    //     let bvh = BVHNode::new(&mut stanford_bunny, 0, size, 0.0, 0.0);

    //     let rays = generate_worse_case_rays(512, &bvh.bounding_box(0.0, 0.0).unwrap());
    //     let mut rays_in = rays.iter().cycle();

    //     bencher.iter(|| {
    //         let ray = rays_in.next().unwrap();
    //         black_box(bvh.hit(&ray, 0.0, f32::INFINITY));
    //     });
    // }

    #[bench]
    fn bench_hit_bvh_stanford_bunny_uniform_rays(bencher: &mut Bencher) {
        let material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
        let mut stanford_bunny = load_hittables_from_obj(Path::new("../input/bunny.obj"), material);
        let size = stanford_bunny.len();
        let bvh = BVHNode::new(&mut stanford_bunny, 0, size, 0.0, 0.0);

        let rays = generate_uniform_random_rays(512, &bvh.bounding_box(0.0, 0.0).unwrap());
        let mut rays_in = rays.iter().cycle();

        bencher.iter(|| {
            let ray = rays_in.next().unwrap();
            black_box(bvh.hit(&ray, 0.0, f32::INFINITY));
        });
    }
}
