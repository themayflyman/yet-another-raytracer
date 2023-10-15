use std::f64::EPSILON;

use crate::aabb::AxisAlignedBoundingBox;
use crate::hittable::{HitRecord, Hittable};
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Clone)]
pub struct Triangle<M: Material> {
    pub vertices: [Vec3; 3],
    pub normals: [Vec3; 3], // array of normal vectors, one per vertex in the mesh. If present, these are interpolated across triangle faces to compute shading normals
    pub uv: [(f64, f64); 3],
    pub material: M,
}

impl<M: Material> Hittable for Triangle<M> {
    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AxisAlignedBoundingBox> {
        let min = Vec3::new(
            self.vertices
                .iter()
                .fold(f64::INFINITY, |min: f64, &v| f64::min(min, v.x())),
            self.vertices
                .iter()
                .fold(f64::INFINITY, |min: f64, &v| f64::min(min, v.y())),
            self.vertices
                .iter()
                .fold(f64::INFINITY, |min: f64, &v| f64::min(min, v.z())),
        );
        let max = Vec3::new(
            self.vertices
                .iter()
                .fold(f64::NEG_INFINITY, |max: f64, &v| f64::max(max, v.x())),
            self.vertices
                .iter()
                .fold(f64::NEG_INFINITY, |max: f64, &v| f64::max(max, v.y())),
            self.vertices
                .iter()
                .fold(f64::NEG_INFINITY, |max: f64, &v| f64::max(max, v.z())),
        );

        return Some(AxisAlignedBoundingBox::new(min, max));
    }

    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
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
