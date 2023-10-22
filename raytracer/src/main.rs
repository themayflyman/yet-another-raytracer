extern crate rand;

use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time;

use indicatif::{ProgressBar, ProgressStyle};

use camera::Camera;
use hittable::{Hittable, HittableList};

use rand::Rng;
use ray::Ray;

use aarect::XZRect;
use color::{gen_wavelength, HasReflectance, CIE_Y_INTERGAL, MAX_LAMBDA, MIN_LAMBDA, RGB, XYZ};
use material::NoMaterial;
use pdf::{HittablePDF, MixurePDF, Pdf};
use scenes::*;
use sphere::StillSphere;
use threadpool::ThreadPool;
use vec3::Vec3;

mod aabb;
mod aarect;
mod box_entity;
mod bvh;
mod camera;
mod color;
mod hittable;
mod material;
mod onb;
mod pdf;
mod ray;
mod scenes;
mod sphere;
mod texture;
mod triangle;
mod vec3;

struct RenderResult {
    pub img: image::RgbaImage,
    pub row: u32,
    pub col: u32,
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

// fn ray_color(
//     r: &Ray,
//     background: RGB,
//     world: &HittableList,
//     lights: Arc<HittableList>,
//     depth: usize,
// ) -> RGB {
//     if depth == 0 {
//         return RGB::default();
//     }
//
//     if let Some(rec) = world.hit(r, 0.001, f64::INFINITY) {
//         let emitted = rec.material.emitted(&rec, rec.u, rec.v, rec.p);
//
//         if let Some(scattered) = rec.material.scatter(r, &rec) {
//             if let Some(scattered_ray) = scattered.ray {
//                 return scattered.color
//                     * ray_color(&scattered_ray, background, world, lights, depth - 1);
//             }
//
//             if scattered.pdf.is_none() {
//                 panic!("Pdf not provided")
//             }
//
//             let mixure_pdf = if lights.objects.is_empty() {
//                 let light_pdf = scattered.pdf.unwrap();
//                 MixurePDF::new(light_pdf.clone(), light_pdf)
//             } else {
//                 let light_pdf = HittablePDF::new(lights.clone(), rec.p);
//                 MixurePDF::new(Arc::new(light_pdf), scattered.pdf.unwrap())
//             };
//             let wl = gen_wavelength();
//             let s = Ray::new(rec.p, mixure_pdf.generate(), r.time(), wl);
//             let pdf_val = mixure_pdf.value(s.direction());
//             return emitted
//                 + scattered.color
//                     * rec.material.scatter_pdf(r, &rec, &s)
//                     * ray_color(&s, background, world, lights, depth - 1)
//                     / pdf_val;
//         } else {
//             return emitted;
//         }
//     }
//
//     background
// }

fn ray_color(
    r: &Ray,
    world: &HittableList,
    lights: Arc<HittableList>,
    background_color: &RGB,
    max_depth: usize,
) -> XYZ {
    let reflectance = ray_reflectance(r, world, lights, background_color, max_depth);
    return XYZ::from_wavelength(r.wavelength) * reflectance;
}

fn ray_reflectance(
    ray_in: &Ray,
    world: &HittableList,
    lights: Arc<HittableList>,
    background_color: &RGB,
    depth: usize,
) -> f64 {
    if depth == 0 {
        return 1.0;
    }

    if let Some(hit_record) = world.hit(ray_in, 0.001, f64::INFINITY) {
        let emitted = hit_record.material.emitted(ray_in, &hit_record);

        if let Some(scattered) = hit_record.material.scatter(ray_in, &hit_record) {
            if let Some(scattered_ray) = scattered.ray {
                return emitted
                    + scattered.attenuation
                        * ray_reflectance(
                            &scattered_ray,
                            world,
                            lights,
                            background_color,
                            depth - 1,
                        );
            }
            if scattered.pdf.is_none() {
                panic!("Pdf not provided")
            }

            let mixure_pdf = if lights.objects.is_empty() {
                let light_pdf = scattered.pdf.unwrap();
                MixurePDF::new(light_pdf.clone(), light_pdf)
            } else {
                let light_pdf = HittablePDF::new(lights.clone(), hit_record.p);
                MixurePDF::new(Arc::new(light_pdf), scattered.pdf.unwrap())
            };
            let wl = gen_wavelength(ray_in.wavelength - 10.0, ray_in.wavelength + 10.0);
            let s = Ray::new(hit_record.p, mixure_pdf.generate(), ray_in.time(), wl);
            let pdf_val = mixure_pdf.value(s.direction(), wl);
            return emitted
                + scattered.attenuation
                    * ray_reflectance(&s, world, lights, background_color, depth - 1)
                    * hit_record.material.scatter_pdf(ray_in, &hit_record, &s)
                    / pdf_val;
        } else {
            return emitted;
        }
    }

    return background_color.reflect(ray_in.wavelength);
}

fn main() {
    // Image
    let mut aspect_ratio: f64 = 3.0 / 2.0;
    let mut image_width: u32 = 1200;
    let mut image_height: u32 = (image_width as f64 / aspect_ratio) as u32;
    let max_depth: usize = 50;

    // World
    let mut world: Arc<HittableList> = Arc::new(HittableList::new());
    let mut lights: HittableList = HittableList::new();
    let lookfrom: Vec3;
    let lookat: Vec3;
    let vfov: f64;
    let mut aperture = 0.0;
    let background: RGB;
    let mut samples_per_pixel: usize = 100;
    let _filename: &str;

    let scene = 11;

    let filename = match scene {
        1 => {
            world = Arc::new(random_scene());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            aperture = 0.1;
            samples_per_pixel = 1000;

            "random_scene.png"
        }

        2 => {
            world = Arc::new(two_spheres());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;

            "two_spheres.png"
        }

        3 => {
            world = Arc::new(two_perlin_spheres());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 30.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;

            "two_perlin_spheres.png"
        }

        4 => {
            world = Arc::new(earth());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;

            "earth.png"
        }

        5 => {
            world = Arc::new(simple_light());
            samples_per_pixel = 400;
            background = RGB::default();
            lookfrom = Vec3::new(26.0, 3.0, 6.0);
            lookat = Vec3::new(0.0, 2.0, 0.0);
            vfov = 20.0;

            "simple_light.png"
        }
        6 => {
            world = Arc::new(cornell_box());
            lights.add_object(Arc::new(XZRect::new(
                213.0, 343.0, 227.0, 332.0, 554.0, NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(190.0, 90.0, 190.0),
                90.0,
                NoMaterial,
            )));
            aspect_ratio = 1.0;
            image_width = 600;
            image_height = 600;
            samples_per_pixel = 100;
            background = RGB::default();
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
            background = RGB::default();
            lookfrom = Vec3::new(278.0, 278.0, -800.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;

            "cornell_box_smoke.png"
        }

        8 => {
            world = Arc::new(the_next_week_final_scene());
            lights.add_object(Arc::new(XZRect::new(
                123.0, 423.0, 147.0, 412.0, 554.0, NoMaterial,
            )));
            aspect_ratio = 1.0;
            image_width = 100;
            image_height = 100;
            samples_per_pixel = 1000;
            background = RGB::default();
            lookfrom = Vec3::new(478.0, 278.0, -600.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;

            "the_next_week_final_scene.png"
        }

        10 => {
            world = Arc::new(simple_light());
            lookfrom = Vec3::new(0.0, 2.0, -10.0);
            lookat = Vec3::new(0.0, 1.0, 0.0);
            aspect_ratio = 2.0;
            image_width = 2000;
            image_height = 1000;
            background = RGB::default();
            samples_per_pixel = 10000;
            aperture = 0.1;
            vfov = 30.0;

            "simple_light.png"
        }

        11 => {
            world = Arc::new(teapot());
            lookfrom = Vec3::new(-10.0, -10.0, 5.0);
            lookat = Vec3::new(0.0, 1.0, 0.0);
            aspect_ratio = 1.0;
            image_width = 200;
            image_height = 200;
            background = RGB::new(0.65, 0.65, 1.0);
            samples_per_pixel = 1000;
            aperture = 0.1;
            vfov = 30.0;

            "teapot.png"
        }

        // 9 => {
        //     world = Arc::new(static_the_next_week_final_scene());
        //     lights.add_object(Arc::new(XZRect::new(
        //         123.0, 423.0, 147.0, 412.0, 554.0, NoMaterial,
        //     )));
        //     aspect_ratio = 1.0;
        //     image_width = 800;
        //     image_height = 800;
        //     samples_per_pixel = 10000;
        //     background = RGB::default();
        //     lookfrom = Vec3::new(478.0, 278.0, -600.0);
        //     lookat = Vec3::new(278.0, 278.0, 0.0);
        //     vfov = 40.0;

        //     "the_next_week_final_static_scene.png"
        // }
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

    let n_workers = 30;
    let mut img = image::RgbaImage::new(image_width as u32, image_height as u32);
    let block_col = 8;
    let block_row = 8;
    let n_jobs = block_col * block_row;
    let pool = ThreadPool::new(n_workers);

    let l = Arc::new(lights);

    let start_time = time::Instant::now();

    let (tx, rx) = channel();
    for col in 0..block_col {
        for row in 0..block_row {
            // let scale: f64 = 1.0 / samples_per_pixel as f64;
            let tx = tx.clone();
            let camera = camera.clone();
            let world = world.clone();
            let lights = l.clone();
            let crop_x = image_width * col / block_col;
            let crop_y = image_height * row / block_row;
            let crop_width = image_width / block_col;
            let crop_height = image_height / block_row;
            let mut imgbuf = image::RgbaImage::new(crop_width, crop_height);

            pool.execute(move || {
                for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
                    // let mut pixel_color = RGB::new(0.0, 0.0, 0.0);
                    // for _s in 0..samples_per_pixel {
                    //     let mut rng = rand::thread_rng();

                    //     let target_x: f64 = x as f64 + crop_x as f64 + rng.gen::<f64>();
                    //     let u: f64 = target_x / (image_width - 1) as f64;
                    //     let target_y: f64 = y as f64 + crop_y as f64 + rng.gen::<f64>();
                    //     let v: f64 = 1.0 - target_y / (image_height - 1) as f64;
                    //     let r: Ray = camera.get_ray(u, v);
                    //     pixel_color = pixel_color
                    //         + ray_color(&r, background, &world, lights.clone(), max_depth);
                    // }

                    // let r = if pixel_color.r().is_nan() {
                    //     0.0
                    // } else {
                    //     pixel_color.r()
                    // };

                    // let g = if pixel_color.g().is_nan() {
                    //     0.0
                    // } else {
                    //     pixel_color.g()
                    // };

                    // let b = if pixel_color.b().is_nan() {
                    //     0.0
                    // } else {
                    //     pixel_color.b()
                    // };

                    // *pixel = image::Rgba([
                    //     (256_f64 * clamp((r * scale).sqrt(), 0.0, 0.999)) as u8,
                    //     (256_f64 * clamp((g * scale).sqrt(), 0.0, 0.999)) as u8,
                    //     (256_f64 * clamp((b * scale).sqrt(), 0.0, 0.999)) as u8,
                    //     255,
                    // ]);
                    let mut pixel_color_xyz = XYZ::default();
                    for _ in 0..samples_per_pixel {
                        let mut rng = rand::thread_rng();
                        let target_x: f64 = x as f64 + crop_x as f64 + rng.gen::<f64>();
                        let u: f64 = target_x / (image_width - 1) as f64;
                        let target_y: f64 = y as f64 + crop_y as f64 + rng.gen::<f64>();
                        let v: f64 = 1.0 - target_y / (image_height - 1) as f64;
                        let wl = gen_wavelength(MIN_LAMBDA, MAX_LAMBDA);
                        let r: Ray = camera.get_ray(u, v, wl);

                        pixel_color_xyz +=
                            ray_color(&r, &world, lights.clone(), &background, max_depth);
                    }

                    pixel_color_xyz = pixel_color_xyz * (MAX_LAMBDA - MIN_LAMBDA)
                        / (CIE_Y_INTERGAL * samples_per_pixel as f64);
                    let pixel_color_rgb = pixel_color_xyz.into_rgb().gamma_corrected();
                    *pixel = image::Rgba([
                        (256_f64 * pixel_color_rgb.r()) as u8,
                        (256_f64 * pixel_color_rgb.g()) as u8,
                        (256_f64 * pixel_color_rgb.b()) as u8,
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
