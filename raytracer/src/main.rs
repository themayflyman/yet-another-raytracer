extern crate rand;

use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time;

use indicatif::{ProgressBar, ProgressStyle};

use camera::Camera;
use hittable::{Hittable, HittableList};

use rand::Rng;
use ray::Ray;

use scenes::*;
use threadpool::ThreadPool;
use vec3::Vec3;

mod aabb;
mod aarect;
mod box_entity;
mod bvh;
mod camera;
mod hittable;
mod material;
mod ray;
mod scenes;
mod sphere;
mod texture;
mod vec3;

struct RenderResult {
    pub img: image::RgbaImage,
    pub row: u32,
    pub col: u32,
}

fn clamp(x: f64, min: f64, max: f64) -> f64 {
    if x < min {
        return min;
    }
    if x > max {
        return max;
    }
    x
}

#[allow(dead_code)]
fn hit_sphere(center: &Vec3, radius: f64, r: &Ray) -> f64 {
    let oc: Vec3 = r.origin() - *center;
    let a: f64 = r.direction().length_squared();
    let half_b: f64 = oc.dot(&r.direction());
    let c: f64 = oc.length_squared() - radius * radius;
    let discriminant: f64 = half_b * half_b - a * c;
    if discriminant < 0.0 {
        -1.0
    } else {
        (-half_b - discriminant.sqrt()) / a
    }
}

fn ray_color(r: &Ray, background: Vec3, world: &HittableList, depth: usize) -> Vec3 {
    if depth == 0 {
        return Vec3::default();
    }

    if let Some(rec) = world.hit(r, 0.001, f64::INFINITY) {
        let scattered = rec.material.scatter(r, &rec);
        let mut emitted = rec.material.emitted(rec.u, rec.v, rec.p);

        if rec.normal.dot(&r.direction()) >= 0.0 {
            emitted = Vec3::default();
        }

        if let Some(scattered_ray) = scattered.ray {
            return emitted
                + scattered.color * ray_color(&scattered_ray, background, world, depth - 1);
        } else {
            return emitted;
        }
    }

    background
}

fn main() {
    // Image
    let mut aspect_ratio: f64 = 16.0 / 9.0;
    let mut image_width: u32 = 400;
    let mut image_height: u32 = (image_width as f64 / aspect_ratio) as u32;
    let max_depth: usize = 50;

    // World
    let world: Arc<HittableList>;
    let lookfrom: Vec3;
    let lookat: Vec3;
    let vfov: f64;
    let mut aperture = 0.0;
    let background: Vec3;
    let mut samples_per_pixel: usize = 100;
    let _filename: &str;

    let scene = 9;

    let filename = match scene {
        1 => {
            world = Arc::new(random_scene());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            aperture = 0.1;

            "random_scene.png"
        }

        2 => {
            world = Arc::new(two_spheres());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;

            "two_spheres.png"
        }

        3 => {
            world = Arc::new(two_perlin_spheres());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;

            "two_perlin_spheres.png"
        }

        4 => {
            world = Arc::new(earth());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;

            "earth.png"
        }

        5 => {
            world = Arc::new(simple_light());
            samples_per_pixel = 400;
            background = Vec3::default();
            lookfrom = Vec3::new(26.0, 3.0, 6.0);
            lookat = Vec3::new(0.0, 2.0, 0.0);
            vfov = 20.0;

            "simple_light.png"
        }

        6 => {
            world = Arc::new(cornell_box());
            aspect_ratio = 1.0;
            image_width = 600;
            image_height = 600;
            samples_per_pixel = 200;
            background = Vec3::default();
            lookfrom = Vec3::new(278.0, 278.0, -800.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;

            "cornell_box.png"
        }

        7 => {
            world = Arc::new(cornell_box_smoke());
            aspect_ratio = 1.0;
            image_width = 600;
            image_height = 600;
            samples_per_pixel = 200;
            background = Vec3::default();
            lookfrom = Vec3::new(278.0, 278.0, -800.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;

            "cornell_box_smoke.png"
        }

        8 => {
            world = Arc::new(the_next_week_final_scene());
            aspect_ratio = 1.0;
            image_width = 800;
            image_height = 800;
            samples_per_pixel = 10000;
            background = Vec3::default();
            lookfrom = Vec3::new(478.0, 278.0, -600.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;

            "the_next_week_final_scene.png"
        }

        9 => {
            world = Arc::new(static_the_next_week_final_scene());
            aspect_ratio = 1.0;
            image_width = 800;
            image_height = 800;
            samples_per_pixel = 10000;
            background = Vec3::default();
            lookfrom = Vec3::new(478.0, 278.0, -600.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;

            "the_next_week_final_static_scene.png"
        }

        _ => {
            panic!("No matched scene");
        }
    };

    let vup: Vec3 = Vec3::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;

    // Camera
    let camera: Arc<Camera> = Arc::new(Camera::new(
        lookfrom,
        lookat,
        vup,
        vfov,
        aspect_ratio,
        aperture,
        dist_to_focus,
        0.0,
        1.0,
    ));

    let n_workers = 8;
    let mut img = image::RgbaImage::new(image_width as u32, image_height as u32);
    let block_col = 8;
    let block_row = 8;
    let n_jobs = block_col * block_row;
    let pool = ThreadPool::new(n_workers);

    let start_time = time::Instant::now();

    let (tx, rx) = channel();
    for col in 0..block_col {
        for row in 0..block_row {
            let scale: f64 = 1.0 / samples_per_pixel as f64;
            let tx = tx.clone();
            let camera = camera.clone();
            let world = world.clone();
            let crop_x = image_width * col / block_col;
            let crop_y = image_height * row / block_row;
            let crop_width = image_width / block_col;
            let crop_height = image_height / block_row;
            let mut imgbuf = image::RgbaImage::new(crop_width, crop_height);

            pool.execute(move || {
                for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
                    let mut pixel_color = Vec3::new(0.0, 0.0, 0.0);
                    for _s in 0..samples_per_pixel {
                        let mut rng = rand::thread_rng();

                        let target_x: f64 = x as f64 + crop_x as f64 + rng.gen::<f64>();
                        let u: f64 = target_x / (image_width - 1) as f64;
                        let target_y: f64 = y as f64 + crop_y as f64 + rng.gen::<f64>();
                        let v: f64 = 1.0 - target_y / (image_height - 1) as f64;
                        let r: Ray = camera.get_ray(u, v);
                        pixel_color = pixel_color + ray_color(&r, background, &world, max_depth);
                    }

                    *pixel = image::Rgba([
                        (256_f64 * clamp((pixel_color.r() * scale).sqrt(), 0.0, 0.999)) as u8,
                        (256_f64 * clamp((pixel_color.g() * scale).sqrt(), 0.0, 0.999)) as u8,
                        (256_f64 * clamp((pixel_color.b() * scale).sqrt(), 0.0, 0.999)) as u8,
                        255,
                    ]);
                }

                tx.send(RenderResult {
                    img: imgbuf,
                    row,
                    col,
                })
                .unwrap();
            })
        }
    }

    let bar = ProgressBar::new(n_jobs as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:>12.cyan.bold} [{bar:57}] {percent}%")
            .progress_chars("=> "),
    );
    bar.set_prefix("Rendering");
    for (_, result) in rx.iter().enumerate().take(n_jobs as usize) {
        let result = result as RenderResult;
        bar.inc(1);

        let x1 = image_width * result.col / block_col;
        let y1 = image_height * result.row / block_row;
        for (x, y, pixel) in result.img.enumerate_pixels() {
            img.put_pixel(x1 + x, y1 + y, *pixel);
        }
    }
    bar.finish_and_clear();

    println!(
        "{} rendered in {} seconds",
        filename,
        start_time.elapsed().as_secs()
    );

    img.save(format!("output/{}", filename)).unwrap();
}
