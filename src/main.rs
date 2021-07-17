mod data;
use data::ray::Ray;
use data::vec3::Vec3;

fn hit_sphere(center: &Vec3, radius: f64, r: &Ray) -> bool {
    let oc: Vec3 = r.origin() - *center;
    let a: f64 = r.direction().dot(&r.direction());
    let b: f64 = 2.0 * oc.dot(&r.direction());
    let c: f64 = oc.dot(&oc) - radius * radius;
    let discriminant: f64 = b * b - 4.0 * a * c;
    discriminant > 0.0
}

// Linear interpolation
fn lerp(r: &Ray) -> Vec3 {
    if hit_sphere(&Vec3::new(0.0, 0.0, -1.0), 0.5, r) {
        return Vec3::new(1.0, 0.0, 0.0);
    }
    let unit_direction: Vec3 = r.direction().unit_vector();
    let t: f64 = 0.5 * (unit_direction.y() + 1.0);
    (1.0 - t) * Vec3::new(1.0, 1.0, 1.0) + t * Vec3::new(0.5, 0.7, 1.0)
}

fn main() {
    // const IMAGE_WIDTH: i32 = 256;
    // const IMAGE_HEIGHT: i32 = 256;

    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: i32 = 400;
    const IMAGE_HEIGHT: i32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as i32;

    const VIEWPORT_HEIGHT: f64 = 2.0;
    const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
    const FOCAL_LENGTH: f64 = 1.0;

    let origin: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    let horizontal: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0.0, 0.0);
    let vertical: Vec3 = Vec3::new(0.0, VIEWPORT_HEIGHT, 0.0);
    let lower_left_corner: Vec3 =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, FOCAL_LENGTH);

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
            let pixel_color: Vec3 = lerp(&r);
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
