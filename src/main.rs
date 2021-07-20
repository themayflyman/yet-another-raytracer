extern crate rand;

use rand::Rng;

use data::ray::Ray;
use data::vec3::Vec3;
use objs::camera::Camera;
use objs::hittable::{HitRecord, Hittable, HittableList};
use objs::sphere::Sphere;

mod data;
mod objs;

fn clamp(x: f64, min: f64, max: f64) -> f64 {
    if x < min {
        return min;
    }
    if x > max {
        return max;
    }
    x
}

fn hit_sphere(center: &Vec3, radius: f64, r: &Ray) -> f64 {
    let oc: Vec3 = r.origin() - *center;
    // let a: f64 = r.direction().dot(&r.direction());
    // let b: f64 = 2.0 * oc.dot(&r.direction());
    // let c: f64 = oc.dot(&oc) - radius * radius;
    // let discriminant: f64 = b * b - 4.0 * a * c;
    // return if discriminant < 0.0 {
    //     -1.0
    // } else {
    //     (-b - discriminant.sqrt()) / (2.0 * a)
    // };
    let a: f64 = r.direction().length_squared();
    let half_b: f64 = oc.dot(&r.direction());
    let c: f64 = oc.length_squared() - radius * radius;
    let discriminant: f64 = half_b * half_b - a * c;
    return if discriminant < 0.0 {
        -1.0
    } else {
        (-half_b - discriminant.sqrt()) / a
    };
}

// Linear interpolation
fn lerp(r: &Ray, world: &HittableList) -> Vec3 {
    let mut rec: HitRecord = HitRecord::new();
    if world.intersect(r, 0.0, f64::INFINITY, &mut rec) {
        return 0.5 * (rec.normal + Vec3::new(1.0, 1.0, 1.0));
    }
    // let mut t: f64 = hit_sphere(&Vec3::new(0.0, 0.0, -1.0), 0.5, r);
    // if t > 0.0 {
    //     let N: Vec3 = (r.at(t) - Vec3::new(0.0, 0.0, -1.0)).unit_vector();
    //     return 0.5 * Vec3::new(N.x() + 1.0, N.y() + 1.0, N.z() + 1.0);
    // }
    let unit_direction: Vec3 = r.direction().unit_vector();
    let t: f64 = 0.5 * (unit_direction.y() + 1.0);
    (1.0 - t) * Vec3::new(1.0, 1.0, 1.0) + t * Vec3::new(0.5, 0.7, 1.0)
}

fn main() {
    // const IMAGE_WIDTH: i32 = 256;
    // const IMAGE_HEIGHT: i32 = 256;

    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: i32 = 400;
    const IMAGE_HEIGHT: i32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as i32;
    const SAMPLES_PER_PIXEL: i32 = 100;

    // World
    let mut world: HittableList = HittableList::new();
    world.add_sphere(Sphere::new(Vec3::new(0.0, 0.0, -1.0), 0.5));
    world.add_sphere(Sphere::new(Vec3::new(0.0, -100.5, -1.0), 100.0));

    // Camera
    let camera: Camera = Camera::new();
    // const VIEWPORT_HEIGHT: f64 = 2.0;
    // const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
    // const FOCAL_LENGTH: f64 = 1.0;
    //
    // let origin: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    // let horizontal: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0.0, 0.0);
    // let vertical: Vec3 = Vec3::new(0.0, VIEWPORT_HEIGHT, 0.0);
    // let lower_left_corner: Vec3 =
    //     origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, FOCAL_LENGTH);

    // Render
    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");

    for j in (0..IMAGE_HEIGHT).rev() {
        eprintln!("Scanlines remaining: {}", j);
        for i in 0..IMAGE_WIDTH {
            let mut pixel_color = Vec3::new(0.0, 0.0, 0.0);
            for _s in 0..SAMPLES_PER_PIXEL {
                let mut rng = rand::thread_rng();

                let u: f64 = (i as f64 + rng.gen::<f64>()) / (IMAGE_WIDTH - 1) as f64;
                let v: f64 = (j as f64 + rng.gen::<f64>()) / (IMAGE_HEIGHT - 1) as f64;
                let r: Ray = camera.get_ray(u, v);
                pixel_color = pixel_color + lerp(&r, &world);
            }

            let scale: f64 = 1.0 / SAMPLES_PER_PIXEL as f64;
            println!(
                "{} {} {}",
                ((256 as f64 * clamp(pixel_color.r() * scale, 0.0, 0.999)) as i32),
                ((256 as f64 * clamp(pixel_color.g() * scale, 0.0, 0.999)) as i32),
                ((256 as f64 * clamp(pixel_color.b() * scale, 0.0, 0.999)) as i32),
            );
        }
    }

    eprintln!("Done.");
}
