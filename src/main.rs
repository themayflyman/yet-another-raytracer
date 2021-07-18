mod data;
mod objs;

use data::ray::Ray;
use data::vec3::Vec3;
use objs::hittable::{HitRecord, Hittable, HittableList};
use objs::sphere::Sphere;

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

    // World
    let mut world: HittableList = HittableList::new();
    world.add_sphere(Sphere::new(Vec3::new(0.0, 0.0, -1.0), 0.5));
    world.add_sphere(Sphere::new(Vec3::new(0.0, -100.5, -1.0), 100.0));

    // Camera
    const VIEWPORT_HEIGHT: f64 = 2.0;
    const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
    const FOCAL_LENGTH: f64 = 1.0;

    let origin: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    let horizontal: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0.0, 0.0);
    let vertical: Vec3 = Vec3::new(0.0, VIEWPORT_HEIGHT, 0.0);
    let lower_left_corner: Vec3 =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, FOCAL_LENGTH);

    // Render
    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");

    for j in (0..IMAGE_HEIGHT).rev() {
        eprintln!("Scanlines remaining: {}", j);
        for i in 0..IMAGE_WIDTH {
            // let r: f64 = i as f64 / (IMAGE_WIDTH - 1) as f64;
            // let g: f64 = j as f64 / (IMAGE_HEIGHT - 1) as f64;
            // let b = 0.25;
            //
            // let ir  = 255.999 * r;
            // let ig = 255.999 * g;
            // let ib = 255.999 * b;
            //
            // let mut color: Vec3 = Vec3::new(ir, ig, ib);
            // color.colorize();

            // println!("{} {} {}", color.r(), color.g(), color.b());
            let u: f64 = i as f64 / (IMAGE_WIDTH - 1) as f64;
            let v: f64 = j as f64 / (IMAGE_HEIGHT - 1) as f64;
            let r: Ray = Ray::new(
                origin,
                lower_left_corner + u * horizontal + v * vertical - origin,
            );
            let pixel_color: Vec3 = lerp(&r, &world);
            println!(
                "{} {} {}",
                (pixel_color.r() * 255.999) as u32,
                (pixel_color.g() * 255.999) as u32,
                (pixel_color.b() * 255.999) as u32
            );
        }
    }

    eprintln!("Done.");
}
