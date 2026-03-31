#![feature(portable_simd)]
#![feature(test)]

extern crate rand;
extern crate test;

use std::path::PathBuf;
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time;

use clap::{Parser, ValueEnum};
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
mod qbvh;
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

enum RenderUpdate {
    Progress(u64),
    Result(RenderResult),
}

const MAX_SAMPLE_LUMINANCE: f64 = 20.0;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum SceneName {
    RandomScene,
    TwoSpheres,
    TwoPerlinSpheres,
    Earth,
    SimpleLight,
    CornellBox,
    CornellBoxSmoke,
    NextWeekFinal,
    Teapot,
    Bunny,
    ThreeSpheres,
    Sycee,
    David,
}

#[derive(Debug, Parser)]
#[command(author, version, about = "Render predefined raytracer scenes")]
struct Cli {
    #[arg(long, value_enum)]
    scene: SceneName,

    #[arg(long)]
    output: Option<PathBuf>,

    #[arg(long, value_parser = clap::value_parser!(u32).range(1..))]
    width: Option<u32>,

    #[arg(long, value_parser = clap::value_parser!(u32).range(1..))]
    height: Option<u32>,

    #[arg(long, value_parser = parse_positive_usize)]
    samples: Option<usize>,

    #[arg(long, value_parser = parse_positive_usize)]
    max_depth: Option<usize>,

    #[arg(long, value_parser = parse_positive_usize)]
    workers: Option<usize>,

    #[arg(long)]
    vfov: Option<f64>,

    #[arg(long)]
    aperture: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct RenderDefaults {
    width: u32,
    height: u32,
    samples_per_pixel: usize,
    max_depth: usize,
    workers: usize,
    vfov: f64,
    aperture: f64,
}

#[derive(Debug, PartialEq)]
struct RenderOptions {
    output_path: PathBuf,
    width: u32,
    height: u32,
    samples_per_pixel: usize,
    max_depth: usize,
    workers: usize,
    vfov: f64,
    aperture: f64,
}

struct ScenePreset {
    defaults: RenderDefaults,
    output_filename: &'static str,
    background: RGB,
    lookfrom: Vec3,
    lookat: Vec3,
    world: Arc<HittableList>,
    lights: HittableList,
}

struct RenderConfig {
    options: RenderOptions,
    background: RGB,
    lookfrom: Vec3,
    lookat: Vec3,
    world: Arc<HittableList>,
    lights: Arc<HittableList>,
}

fn parse_positive_usize(input: &str) -> Result<usize, String> {
    let value = input
        .parse::<usize>()
        .map_err(|err| format!("invalid integer `{input}`: {err}"))?;
    if value == 0 {
        return Err("value must be greater than 0".to_string());
    }

    Ok(value)
}

fn default_output_path(filename: &str) -> PathBuf {
    PathBuf::from("output").join(filename)
}

fn resolve_dimensions(
    default_width: u32,
    default_height: u32,
    width_override: Option<u32>,
    height_override: Option<u32>,
) -> (u32, u32) {
    let aspect_ratio = default_width as f64 / default_height as f64;

    match (width_override, height_override) {
        (Some(width), Some(height)) => (width, height),
        (Some(width), None) => {
            let height = ((width as f64 / aspect_ratio).round().max(1.0)) as u32;
            (width, height)
        }
        (None, Some(height)) => {
            let width = ((height as f64 * aspect_ratio).round().max(1.0)) as u32;
            (width, height)
        }
        (None, None) => (default_width, default_height),
    }
}

fn resolve_render_options(
    default_filename: &str,
    defaults: RenderDefaults,
    cli: &Cli,
) -> RenderOptions {
    let (width, height) =
        resolve_dimensions(defaults.width, defaults.height, cli.width, cli.height);

    RenderOptions {
        output_path: cli
            .output
            .clone()
            .unwrap_or_else(|| default_output_path(default_filename)),
        width,
        height,
        samples_per_pixel: cli.samples.unwrap_or(defaults.samples_per_pixel),
        max_depth: cli.max_depth.unwrap_or(defaults.max_depth),
        workers: cli.workers.unwrap_or(defaults.workers),
        vfov: cli.vfov.unwrap_or(defaults.vfov),
        aperture: cli.aperture.unwrap_or(defaults.aperture),
    }
}

fn build_scene_preset(scene: SceneName) -> ScenePreset {
    let mut defaults = RenderDefaults {
        width: 1200,
        height: 800,
        samples_per_pixel: 100,
        max_depth: 50,
        workers: 30,
        vfov: 20.0,
        aperture: 0.0,
    };
    let mut lights = HittableList::new();
    let world: Arc<HittableList>;
    let lookfrom: Vec3;
    let lookat: Vec3;
    let background: RGB;
    let output_filename: &'static str;

    match scene {
        SceneName::RandomScene => {
            world = Arc::new(random_scene());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            defaults.aperture = 0.1;
            defaults.samples_per_pixel = 1000;
            output_filename = "random_scene.png";
        }
        SceneName::TwoSpheres => {
            world = Arc::new(two_spheres());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            output_filename = "two_spheres.png";
        }
        SceneName::TwoPerlinSpheres => {
            world = Arc::new(two_perlin_spheres());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 30.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            output_filename = "two_perlin_spheres.png";
        }
        SceneName::Earth => {
            world = Arc::new(earth());
            background = RGB::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            output_filename = "earth.png";
        }
        SceneName::SimpleLight => {
            world = Arc::new(simple_light());
            background = RGB::default();
            lookfrom = Vec3::new(26.0, 3.0, 6.0);
            lookat = Vec3::new(0.0, 2.0, 0.0);
            defaults.samples_per_pixel = 400;
            output_filename = "simple_light.png";
        }
        SceneName::CornellBox => {
            world = Arc::new(cornell_box());
            lights.add_object(Arc::new(XZRect::new(
                213.0, 343.0, 227.0, 332.0, 554.0, NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(190.0, 90.0, 190.0),
                90.0,
                NoMaterial,
            )));
            background = RGB::default();
            lookfrom = Vec3::new(278.0, 278.0, -800.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            defaults.width = 600;
            defaults.height = 600;
            defaults.samples_per_pixel = 100;
            defaults.vfov = 40.0;
            output_filename = "cornell_box.png";
        }
        SceneName::CornellBoxSmoke => {
            world = Arc::new(cornell_box_smoke());
            background = RGB::default();
            lookfrom = Vec3::new(278.0, 278.0, -800.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            defaults.width = 600;
            defaults.height = 600;
            defaults.samples_per_pixel = 200;
            defaults.vfov = 40.0;
            output_filename = "cornell_box_smoke.png";
        }
        SceneName::NextWeekFinal => {
            world = Arc::new(the_next_week_final_scene());
            lights.add_object(Arc::new(XZRect::new(
                123.0, 423.0, 147.0, 412.0, 554.0, NoMaterial,
            )));
            background = RGB::default();
            lookfrom = Vec3::new(478.0, 278.0, -600.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            defaults.width = 100;
            defaults.height = 100;
            defaults.samples_per_pixel = 1000;
            defaults.vfov = 40.0;
            output_filename = "the_next_week_final_scene.png";
        }
        SceneName::Teapot => {
            world = Arc::new(teapot());
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(30.0, 40.0, -30.0),
                20.0,
                NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(-20.0, 10.0, 50.0),
                10.0,
                NoMaterial,
            )));
            background = RGB::new(0.0, 0.0, 0.0);
            lookfrom = Vec3::new(5.0, 50.0, 60.0);
            lookat = Vec3::new(0.0, 5.0, 0.0);
            defaults.width = 1000;
            defaults.height = 1000;
            defaults.samples_per_pixel = 8000;
            defaults.vfov = 30.0;
            defaults.aperture = 0.001;
            output_filename = "teapot.png";
        }
        SceneName::Bunny => {
            world = Arc::new(bunny());
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(0.0, 6.0, 2.0),
                2.0,
                NoMaterial,
            )));
            background = RGB::new(0.0, 0.0, 0.0);
            lookfrom = Vec3::new(0.0, 2.0, 10.0);
            lookat = Vec3::new(0.0, 1.0, 0.0);
            defaults.width = 1000;
            defaults.height = 1000;
            defaults.samples_per_pixel = 50;
            defaults.vfov = 30.0;
            defaults.aperture = 0.1;
            output_filename = "bunny.png";
        }
        SceneName::ThreeSpheres => {
            world = Arc::new(three_spheres());
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(0.0, 6.0, 2.0),
                2.0,
                NoMaterial,
            )));
            background = RGB::new(0.0, 0.0, 0.0);
            lookfrom = Vec3::new(1.0, 5.0, -8.0);
            lookat = Vec3::new(0.0, 1.0, 0.0);
            defaults.width = 1000;
            defaults.height = 1000;
            defaults.samples_per_pixel = 5000;
            defaults.vfov = 30.0;
            defaults.aperture = 0.1;
            output_filename = "three_spheres.png";
        }
        SceneName::Sycee => {
            world = Arc::new(sycee());
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(0.0, 6.0, 2.0),
                2.0,
                NoMaterial,
            )));
            background = RGB::new(0.0, 0.0, 0.0);
            lookfrom = Vec3::new(1.0, 5.0, -8.0);
            lookat = Vec3::new(0.0, 1.0, 0.0);
            defaults.width = 1000;
            defaults.height = 1000;
            defaults.samples_per_pixel = 5000;
            defaults.vfov = 30.0;
            defaults.aperture = 0.1;
            output_filename = "sycee.png";
        }
        SceneName::David => {
            world = Arc::new(david());
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(1200.0, 1300.0, 800.0),
                700.0,
                NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(-1200.0, 1300.0, 800.0),
                700.0,
                NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(1200.0, 1300.0, -800.0),
                700.0,
                NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(1200.0, -1300.0, -800.0),
                700.0,
                NoMaterial,
            )));
            lights.add_object(Arc::new(StillSphere::new(
                Vec3::new(1200.0, 1300.0, -800.0),
                700.0,
                NoMaterial,
            )));
            background = RGB::new(0.0, 0.0, 0.0);
            lookfrom = Vec3::new(50.0, 120.0, 300.0);
            lookat = Vec3::new(0.0, 120.0, 0.0);
            defaults.width = 600;
            defaults.height = 600;
            defaults.samples_per_pixel = 10000;
            defaults.vfov = 20.0;
            defaults.aperture = 0.001;
            output_filename = "david.png";
        }
    }

    ScenePreset {
        defaults,
        output_filename,
        background,
        lookfrom,
        lookat,
        world,
        lights,
    }
}

fn resolve_render_config(cli: Cli) -> RenderConfig {
    let preset = build_scene_preset(cli.scene);
    let options = resolve_render_options(preset.output_filename, preset.defaults, &cli);

    RenderConfig {
        options,
        background: preset.background,
        lookfrom: preset.lookfrom,
        lookat: preset.lookat,
        world: preset.world,
        lights: Arc::new(preset.lights),
    }
}

fn sanitize_sample_xyz(sample: XYZ) -> XYZ {
    if !sample.x().is_finite() || !sample.y().is_finite() || !sample.z().is_finite() {
        return XYZ::default();
    }

    let luminance = sample.y();
    if luminance <= 0.0 || luminance <= MAX_SAMPLE_LUMINANCE {
        return sample;
    }

    sample * (MAX_SAMPLE_LUMINANCE / luminance)
}

fn clamp_display_channel(channel: f64) -> u8 {
    (256.0 * channel.clamp(0.0, 0.999)) as u8
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
                return scattered.attenuation
                    * ray_reflectance(&scattered_ray, world, lights, background_color, depth - 1);
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
            let s = Ray::new(
                hit_record.p,
                mixure_pdf.generate(),
                ray_in.time(),
                ray_in.wavelength,
            );
            let pdf_val = mixure_pdf.value(s.direction(), ray_in.wavelength);
            if !pdf_val.is_finite() || pdf_val <= 0.0 {
                return emitted;
            }

            return scattered.attenuation
                * ray_reflectance(&s, world, lights, background_color, depth - 1)
                * hit_record.material.scatter_pdf(ray_in, &hit_record, &s)
                / pdf_val;
        } else {
            return emitted;
        }
    }

    return background_color.reflect(ray_in.wavelength);
}

fn render(config: RenderConfig) {
    let start_time = time::Instant::now();
    let RenderConfig {
        options,
        background,
        lookfrom,
        lookat,
        world,
        lights,
    } = config;
    let RenderOptions {
        output_path,
        width: image_width,
        height: image_height,
        samples_per_pixel,
        max_depth,
        workers,
        vfov,
        aperture,
    } = options;
    let aspect_ratio = image_width as f64 / image_height as f64;

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

    let mut img = image::RgbaImage::new(image_width as u32, image_height as u32);
    let block_col = 8;
    let block_row = 8;
    let n_jobs = block_col * block_row;
    let total_pixels = u64::from(image_width) * u64::from(image_height);
    let pool = ThreadPool::new(workers);

    let (tx, rx) = channel::<RenderUpdate>();
    for col in 0..block_col {
        for row in 0..block_row {
            // let scale: f64 = 1.0 / samples_per_pixel as f64;
            let tx = tx.clone();
            let camera = camera.clone();
            let world = world.clone();
            let lights = lights.clone();
            let crop_x = image_width * col / block_col;
            let crop_y = image_height * row / block_row;
            let crop_width = image_width / block_col;
            let crop_height = image_height / block_row;
            let mut imgbuf = image::RgbaImage::new(crop_width, crop_height);

            pool.execute(move || {
                for y in 0..crop_height {
                    for x in 0..crop_width {
                        let pixel = imgbuf.get_pixel_mut(x, y);
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

                            let sample_xyz = sanitize_sample_xyz(ray_color(
                                &r,
                                &world,
                                lights.clone(),
                                &background,
                                max_depth,
                            ));
                            pixel_color_xyz += sample_xyz;
                        }

                        pixel_color_xyz = pixel_color_xyz * (MAX_LAMBDA - MIN_LAMBDA)
                            / (CIE_Y_INTERGAL * samples_per_pixel as f64);
                        let pixel_color_rgb = pixel_color_xyz.into_rgb().gamma_corrected();
                        *pixel = image::Rgba([
                            clamp_display_channel(pixel_color_rgb.r()),
                            clamp_display_channel(pixel_color_rgb.g()),
                            clamp_display_channel(pixel_color_rgb.b()),
                            255,
                        ]);
                    }
                    tx.send(RenderUpdate::Progress(u64::from(crop_width)))
                        .unwrap();
                }

                tx.send(RenderUpdate::Result(RenderResult {
                    img: imgbuf,
                    row,
                    col,
                }))
                .unwrap();
            })
        }
    }
    drop(tx);

    let bar = ProgressBar::new(total_pixels);
    bar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{prefix:>12.cyan.bold} [{elapsed_precise}] [{bar:28.cyan/blue}] \
                 {percent:>3}% {pos:>7}/{len:7} px {per_sec:>8} ETA {eta_precise}",
            )
            .expect("progress bar template should be valid")
            .progress_chars("=> "),
    );
    bar.set_prefix("Rendering");
    let mut completed_jobs = 0;
    while completed_jobs < n_jobs as usize {
        match rx.recv().unwrap() {
            RenderUpdate::Progress(pixels) => bar.inc(pixels),
            RenderUpdate::Result(result) => {
                completed_jobs += 1;

                let x1 = image_width * result.col / block_col;
                let y1 = image_height * result.row / block_row;
                for (x, y, pixel) in result.img.enumerate_pixels() {
                    img.put_pixel(x1 + x, y1 + y, *pixel);
                }
            }
        }
    }
    bar.finish_and_clear();

    println!(
        "{} rendered in {} seconds",
        output_path.display(),
        start_time.elapsed().as_secs()
    );

    if let Some(parent) = output_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).unwrap();
        }
    }
    img.save(&output_path).unwrap();
}

fn main() {
    let cli = Cli::parse();
    let config = resolve_render_config(cli);
    render(config);
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use clap::Parser;

    use super::{
        default_output_path, resolve_dimensions, resolve_render_options, sanitize_sample_xyz, Cli,
        RenderDefaults, SceneName, MAX_SAMPLE_LUMINANCE, XYZ,
    };

    fn test_cli() -> Cli {
        Cli {
            scene: SceneName::David,
            output: None,
            width: None,
            height: None,
            samples: None,
            max_depth: None,
            workers: None,
            vfov: None,
            aperture: None,
        }
    }

    #[test]
    fn sanitize_sample_xyz_drops_non_finite_samples() {
        assert_eq!(
            sanitize_sample_xyz(XYZ::new(f64::NAN, 1.0, 1.0)),
            XYZ::default()
        );
        assert_eq!(
            sanitize_sample_xyz(XYZ::new(f64::INFINITY, 1.0, 1.0)),
            XYZ::default()
        );
    }

    #[test]
    fn sanitize_sample_xyz_clamps_luminance_but_preserves_chromaticity() {
        let sample = XYZ::new(40.0, 80.0, 20.0);
        let clamped = sanitize_sample_xyz(sample);

        assert!((clamped.y() - MAX_SAMPLE_LUMINANCE).abs() < 1e-9);
        assert!((clamped.x() / clamped.y() - sample.x() / sample.y()).abs() < 1e-9);
        assert!((clamped.z() / clamped.y() - sample.z() / sample.y()).abs() < 1e-9);
    }

    #[test]
    fn cli_accepts_named_scene_values() {
        let cli = Cli::try_parse_from(["raytracer", "--scene", "david"]).unwrap();

        assert_eq!(cli.scene, SceneName::David);
    }

    #[test]
    fn cli_rejects_unknown_scene_values() {
        let err = Cli::try_parse_from(["raytracer", "--scene", "unknown-scene"]).unwrap_err();

        assert_eq!(err.kind(), clap::error::ErrorKind::InvalidValue);
    }

    #[test]
    fn resolve_dimensions_uses_defaults_without_overrides() {
        assert_eq!(resolve_dimensions(1200, 800, None, None), (1200, 800));
    }

    #[test]
    fn resolve_dimensions_derives_height_from_width() {
        assert_eq!(resolve_dimensions(1200, 800, Some(600), None), (600, 400));
    }

    #[test]
    fn resolve_dimensions_derives_width_from_height() {
        assert_eq!(resolve_dimensions(1200, 800, None, Some(400)), (600, 400));
    }

    #[test]
    fn resolve_dimensions_uses_explicit_width_and_height() {
        assert_eq!(
            resolve_dimensions(1200, 800, Some(1024), Some(512)),
            (1024, 512)
        );
    }

    #[test]
    fn resolve_render_options_uses_default_output_when_not_overridden() {
        let defaults = RenderDefaults {
            width: 1200,
            height: 800,
            samples_per_pixel: 100,
            max_depth: 50,
            workers: 30,
            vfov: 20.0,
            aperture: 0.0,
        };
        let cli = test_cli();

        let options = resolve_render_options("david.png", defaults, &cli);

        assert_eq!(options.output_path, default_output_path("david.png"));
    }

    #[test]
    fn resolve_render_options_respects_output_and_scalar_overrides() {
        let defaults = RenderDefaults {
            width: 1200,
            height: 800,
            samples_per_pixel: 100,
            max_depth: 50,
            workers: 30,
            vfov: 20.0,
            aperture: 0.0,
        };
        let mut cli = test_cli();
        cli.output = Some(PathBuf::from("custom/output.png"));
        cli.width = Some(600);
        cli.samples = Some(32);
        cli.max_depth = Some(12);
        cli.workers = Some(8);
        cli.vfov = Some(45.0);
        cli.aperture = Some(0.25);

        let options = resolve_render_options("david.png", defaults, &cli);

        assert_eq!(options.output_path, PathBuf::from("custom/output.png"));
        assert_eq!(options.width, 600);
        assert_eq!(options.height, 400);
        assert_eq!(options.samples_per_pixel, 32);
        assert_eq!(options.max_depth, 12);
        assert_eq!(options.workers, 8);
        assert_eq!(options.vfov, 45.0);
        assert_eq!(options.aperture, 0.25);
    }
}
