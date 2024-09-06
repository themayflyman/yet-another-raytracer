use std::f32::EPSILON;
use std::path::Path;
use std::sync::Arc;

use crate::aabb::AxisAlignedBoundingBox;
use crate::bvh::BVHNode;
use crate::hittable::{HitRecord, Hittable, HittableList};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct Triangle<M: Material> {
    pub vertices: [Vec3; 3],
    pub normals: [Vec3; 3], // array of normal vectors, one per vertex in the mesh. If present, these are interpolated across triangle faces to compute shading normals
    pub uv: [(f32, f32); 3],
    pub material: M,
}

impl<M: Material> Hittable for Triangle<M> {
    fn bounding_box(&self, _time0: f32, _time1: f32) -> Option<AxisAlignedBoundingBox> {
        let min = Vec3::new(
            self.vertices
                .iter()
                .fold(f32::INFINITY, |min: f32, &v| f32::min(min, v.x())),
            self.vertices
                .iter()
                .fold(f32::INFINITY, |min: f32, &v| f32::min(min, v.y())),
            self.vertices
                .iter()
                .fold(f32::INFINITY, |min: f32, &v| f32::min(min, v.z())),
        );
        let max = Vec3::new(
            self.vertices
                .iter()
                .fold(f32::NEG_INFINITY, |max: f32, &v| f32::max(max, v.x())),
            self.vertices
                .iter()
                .fold(f32::NEG_INFINITY, |max: f32, &v| f32::max(max, v.y())),
            self.vertices
                .iter()
                .fold(f32::NEG_INFINITY, |max: f32, &v| f32::max(max, v.z())),
        );

        return Some(AxisAlignedBoundingBox::new(min, max));
    }

    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let edge1 = self.vertices[1] - self.vertices[0];
        let edge2 = self.vertices[2] - self.vertices[0];

        let h = ray.direction().cross(&edge2);
        let a = edge1.dot(&h);

        if a > -EPSILON && a < EPSILON {
            return None; // This ray is parallel to this triangle.
        }

        let f = 1.0 / a;
        let s = ray.origin() - self.vertices[0];
        let u = f * s.dot(&h);

        if u < 0.0 || u > 1.0 {
            return None;
        }

        let q = s.cross(&edge1);
        let v = f * ray.direction().dot(&q);

        if v < 0.0 || u + v > 1.0 {
            return None;
        }

        // At this stage we can compute t to find out where the intersection point is on the line.
        let t = f * edge2.dot(&q);
        if t < t_min || t > t_max {
            // This means that there is a line intersection but not a ray intersection.
            return None;
        }
        let p: Vec3 = ray.at(t);
        let w = 1.0 - u - v;
        let outward_normal =
            self.normals[0] * v + self.normals[1] * u + self.normals[2] * w;
        let normal: Vec3;
        let front_face: bool;
        if ray.direction().dot(&outward_normal) < 0.0 {
            normal = outward_normal;
            front_face = true;
        } else {
            normal = -outward_normal;
            front_face = false;
        }
        return Some(HitRecord {
            u: self.uv[0].0 * v + self.uv[1].0 * u + self.uv[2].0 * w,
            v: self.uv[0].1 * v + self.uv[1].1 * u + self.uv[2].1 * w,
            t,
            p,
            normal,
            front_face,
            material: &self.material,
        });
    }
}

#[derive(Clone)]
pub struct TriangleMesh {
    // to accelerate triangle hit test
    pub bvh: BVHNode,
}

impl TriangleMesh {
    pub fn from_obj(obj_path: &Path, material: impl Material + Clone + 'static) -> Self {
        let (models, _materials) =
            tobj::load_obj(obj_path, &tobj::GPU_LOAD_OPTIONS).expect("Failed to load OBJ file");

        let mut read_mesh = HittableList::new();
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
                let default_normal = (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0])).unit_vector();
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
                read_mesh.add_object(Arc::new(triangle));
            }
        }

        let size = read_mesh.size();
        let bvh = BVHNode::new(&mut read_mesh.objects, 0, size, 0.0, 0.0);

        return Self { bvh };
    }
}

impl Hittable for TriangleMesh {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        return self.bvh.hit(ray, t_min, t_max);
    }

    fn bounding_box(&self, time0: f32, time1: f32) -> Option<AxisAlignedBoundingBox> {
        return self.bvh.bounding_box(time0, time1);
    }
}
