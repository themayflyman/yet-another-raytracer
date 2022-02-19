extern crate rand;

use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time;

use aarect::{XYRect, XZRect, YZRect};
use bvh::BVHNode;
use camera::Camera;
use hittable::{Hittable, HittableList};
use material::{Dielectric, DiffuseLight, Lambertian, Metal};
use rand::Rng;
use ray::Ray;
use sphere::{MovingSphere, StillSphere};
use texture::{CheckerTexture, NoiseTexture, SolidColor};
use texture::{ImageTexture, NoiseType};
use threadpool::ThreadPool;
use vec3::Vec3;

use self::box_entity::BoxEntity;
use self::hittable::{ConstantMedium, RotateY, Translate};

mod aabb;
mod aarect;
mod box_entity;
mod bvh;
mod camera;
mod hittable;
mod material;
mod ray;
mod sphere;
mod texture;
mod vec3;

struct RenderResult {
    pub pixel: image::Rgba<u8>,
    pub x: u32,
    pub y: u32,
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
    return if discriminant < 0.0 {
        -1.0
    } else {
        (-half_b - discriminant.sqrt()) / a
    };
}

fn ray_color(r: &Ray, background: Vec3, world: &HittableList, depth: usize) -> Vec3 {
    if depth <= 0 {
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

    return background;
}

fn random_scene() -> HittableList {
    let mut world: HittableList = HittableList::new();

    let checker = CheckerTexture::new(Vec3::new(0.2, 0.3, 0.1), Vec3::new(0.9, 0.9, 0.9));
    let group_material = Lambertian::new(checker);
    world.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        group_material,
    )));

    for a in -11..11 {
        for b in -11..11 {
            let mut rng = rand::thread_rng();
            let choose_mat: f64 = rng.gen::<f64>();
            let center: Vec3 = Vec3::new(
                a as f64 + 0.9 * rng.gen::<f64>(),
                0.2,
                b as f64 + 0.9 * rng.gen::<f64>(),
            );

            if (center - Vec3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                if choose_mat < 0.8 {
                    let albedo = SolidColor::new(Vec3::new(
                        rng.gen_range(-1.0, 1.0),
                        rng.gen_range(-1.0, 1.0),
                        rng.gen_range(-1.0, 1.0),
                    ));
                    let sphere_material = Lambertian::new(albedo);
                    let center2 = center + Vec3::new(0.0, rng.gen_range(0.0, 0.5), 0.0);
                    world.add_sphere(Arc::new(MovingSphere::new(
                        center,
                        center2,
                        0.0,
                        1.0,
                        0.2,
                        sphere_material,
                    )));
                } else if choose_mat < 0.95 {
                    let albedo: Vec3 = Vec3::new(
                        rng.gen_range(0.5, 1.0),
                        rng.gen_range(0.5, 1.0),
                        rng.gen_range(0.5, 1.0),
                    );
                    let fuzz: f64 = rng.gen_range(0.0, 0.5);
                    let sphere_material = Metal::new(albedo, fuzz);
                    world.add_sphere(Arc::new(StillSphere::new(center, 0.2, sphere_material)));
                } else {
                    let sphere_material = Dielectric::new(1.5);
                    world.add_sphere(Arc::new(StillSphere::new(center, 0.2, sphere_material)));
                }
            }
        }
    }

    world.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, 1.0, 0.0),
        1.0,
        Dielectric::new(1.5),
    )));

    world.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(-4.0, 1.0, 0.0),
        1.0,
        Dielectric::new(1.5),
    )));

    world.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(4.0, 1.0, 0.0),
        1.0,
        Dielectric::new(1.5),
    )));

    return world;
}

fn two_spheres() -> HittableList {
    let mut objects = HittableList::new();

    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, -10.0, 0.0),
        10.0,
        Lambertian::new(CheckerTexture::new(
            Vec3::new(0.2, 0.3, 0.1),
            Vec3::new(0.9, 0.9, 0.9),
        )),
    )));
    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, 10.0, 0.0),
        10.0,
        Lambertian::new(CheckerTexture::new(
            Vec3::new(0.2, 0.3, 0.1),
            Vec3::new(0.9, 0.9, 0.9),
        )),
    )));

    objects
}

fn two_perlin_spheres() -> HittableList {
    let mut objects = HittableList::new();

    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));
    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, 2.0, 0.0),
        2.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));

    objects
}

fn earth() -> HittableList {
    let mut objects = HittableList::new();

    let earth_texture = ImageTexture::new("earthmap.jpg").unwrap();
    let earth_surface = Lambertian::new(earth_texture);
    let global = StillSphere::new(Vec3::new(0.0, 0.0, 0.0), 2.0, earth_surface);
    objects.add_sphere(Arc::new(global));

    objects
}

fn simple_light() -> HittableList {
    let mut objects = HittableList::new();

    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));
    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, 2.0, 0.0),
        2.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));

    let difflight = DiffuseLight::new(SolidColor::new(Vec3::new(4.0, 4.0, 4.0)));
    objects.add_sphere(Arc::new(XYRect::new(3.0, 5.0, 1.0, 3.0, -2.0, difflight)));

    objects
}

fn cornell_box() -> HittableList {
    let mut objects = HittableList::new();

    let red = Lambertian::new(SolidColor::new(Vec3::new(0.65, 0.05, 0.05)));
    let white = Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73)));
    let green = Lambertian::new(SolidColor::new(Vec3::new(0.12, 0.45, 0.15)));
    let light = DiffuseLight::new(SolidColor::new(Vec3::new(15.0, 15.0, 15.0)));

    objects.add_sphere(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 555.0, green)));
    objects.add_sphere(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 0.0, red)));
    objects.add_sphere(Arc::new(XZRect::new(
        213.0, 343.0, 227.0, 332.0, 554.0, light,
    )));
    objects.add_sphere(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73))),
    )));
    objects.add_sphere(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73))),
    )));
    objects.add_sphere(Arc::new(XYRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73))),
    )));

    let box1 = Arc::new(Translate::new(
        RotateY::new(
            BoxEntity::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(165.0, 330.0, 165.0),
                Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73))),
            ),
            15.0,
        ),
        Vec3::new(265.0, 0.0, 295.0),
    ));
    let box2 = Arc::new(Translate::new(
        RotateY::new(
            BoxEntity::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(165.0, 165.0, 165.0),
                Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73))),
            ),
            -18.0,
        ),
        Vec3::new(130.0, 0.0, 65.0),
    ));

    objects.add_sphere(box1);
    objects.add_sphere(box2);

    objects
}

fn cornell_box_smoke() -> HittableList {
    let mut objects = HittableList::new();

    let red = Lambertian::new(SolidColor::new(Vec3::new(0.65, 0.05, 0.05)));
    let white = Lambertian::new(SolidColor::new(Vec3::new(0.73, 0.73, 0.73)));
    let green = Lambertian::new(SolidColor::new(Vec3::new(0.12, 0.45, 0.15)));
    let light = DiffuseLight::new(SolidColor::new(Vec3::new(7.0, 7.0, 7.0)));

    objects.add_sphere(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 555.0, green)));
    objects.add_sphere(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 0.0, red)));
    objects.add_sphere(Arc::new(XZRect::new(
        113.0, 443.0, 127.0, 432.0, 554.0, light,
    )));
    objects.add_sphere(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        white.clone(),
    )));
    objects.add_sphere(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));
    objects.add_sphere(Arc::new(XYRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));

    let box1 = Arc::new(ConstantMedium::new(
        Translate::new(
            RotateY::new(
                BoxEntity::new(
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(165.0, 330.0, 165.0),
                    white.clone(),
                ),
                15.0,
            ),
            Vec3::new(265.0, 0.0, 295.0),
        ),
        0.01,
        SolidColor::new(Vec3::new(0.0, 0.0, 0.0)),
    ));
    let box2 = Arc::new(ConstantMedium::new(
        Translate::new(
            RotateY::new(
                BoxEntity::new(
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(165.0, 165.0, 165.0),
                    white.clone(),
                ),
                -18.0,
            ),
            Vec3::new(130.0, 0.0, 65.0),
        ),
        0.01,
        SolidColor::new(Vec3::new(1.0, 1.0, 1.0)),
    ));

    objects.add_sphere(box1);
    objects.add_sphere(box2);

    objects
}

fn the_next_week_final_scene() -> HittableList {
    let mut rng = rand::thread_rng();

    let mut boxes1 = HittableList::new();
    let ground = Lambertian::new(SolidColor::new(Vec3::new(0.48, 0.83, 0.53)));

    for i in 0..20 {
        for j in 0..20 {
            let w = 100.0;
            let x0 = -1000.0 + i as f64 * w;
            let z0 = -1000.0 + j as f64 * w;
            let y0 = 0.0;
            let x1 = x0 + w;
            let y1 = rng.gen_range(1.0, 101.0);
            let z1 = z0 + w;

            boxes1.add_sphere(Arc::new(BoxEntity::new(
                Vec3::new(x0, y0, z0),
                Vec3::new(x1, y1, z1),
                ground.clone(),
            )));
        }
    }

    let mut objects = HittableList::new();

    objects.add_sphere(Arc::new(BVHNode::new(&mut boxes1.spheres, 0, 400, 0.0, 0.0)));

    let light = DiffuseLight::new(SolidColor::new(Vec3::new(7.0, 7.0, 7.0)));
    objects.add_sphere(Arc::new(XZRect::new(
        123.0, 423.0, 147.0, 412.0, 554.0, light,
    )));

    let center0 = Vec3::new(400.0, 400.0, 200.0);
    let center1 = center0 + Vec3::new(30.0, 0.0, 0.0);
    let moving_sphere_material = Lambertian::new(SolidColor::new(Vec3::new(0.7, 0.3, 0.1)));
    objects.add_sphere(Arc::new(MovingSphere::new(
        center0,
        center1,
        0.0,
        1.0,
        50.0,
        moving_sphere_material,
    )));

    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(260.0, 150.0, 45.0),
        50.0,
        Dielectric::new(1.5),
    )));
    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(0.0, 150.0, 145.0),
        50.0,
        Metal::new(Vec3::new(0.8, 0.8, 0.9), 1.0),
    )));

    let mut boundary = StillSphere::new(Vec3::new(360.0, 150.0, 145.0), 70.0, Dielectric::new(1.5));
    let boundary2 = StillSphere::new(Vec3::new(360.0, 150.0, 145.0), 70.0, Dielectric::new(1.5));
    objects.add_sphere(Arc::new(boundary));
    objects.add_sphere(Arc::new(ConstantMedium::new(
        boundary2,
        0.2,
        SolidColor::new(Vec3::new(0.2, 0.4, 0.9)),
    )));
    boundary = StillSphere::new(Vec3::new(0.0, 0.0, 0.0), 5000.0, Dielectric::new(1.5));
    objects.add_sphere(Arc::new(ConstantMedium::new(
        boundary,
        0.0001,
        SolidColor::new(Vec3::new(1.0, 1.0, 1.0)),
    )));

    let emat = Lambertian::new(ImageTexture::new("earthmap.jpg").unwrap());
    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(400.0, 200.0, 400.0),
        100.0,
        emat,
    )));
    let pertext = NoiseTexture::new(NoiseType::Marble, 0.1);
    objects.add_sphere(Arc::new(StillSphere::new(
        Vec3::new(220.0, 280.0, 300.0),
        80.0,
        Lambertian::new(pertext),
    )));

    let mut boxes2 = HittableList::new();
    let white = Lambertian::new(SolidColor::new_from_value(0.73, 0.73, 0.73));
    for j in 0..1000 {
        boxes2.add_sphere(Arc::new(StillSphere::new(
            Vec3::random(0.0, 165.0),
            10.0,
            white.clone(),
        )));
    }

    objects.add_sphere(Arc::new(Translate::new(
        RotateY::new(BVHNode::new(&mut boxes2.spheres, 0, 1000, 0.0, 0.0), 15.0),
        Vec3::new(-100.0, 270.0, 395.0),
    )));

    objects
}

fn main() {
    // Image
    let mut aspect_ratio: f64 = 16.0 / 9.0;
    let mut image_width: u32 = 400;
    let mut image_height: u32 = (image_width as f64 / aspect_ratio) as u32;
    let mut max_depth: usize = 50;

    // World
    let world: Arc<HittableList>;
    let lookfrom: Vec3;
    let lookat: Vec3;
    let vfov: f64;
    let mut aperture = 0.0;
    let background: Vec3;
    let mut samples_per_pixel: usize = 100;
    let mut filename = "default.png";

    let scene = 8;

    match scene {
        1 => {
            world = Arc::new(random_scene());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            aperture = 0.1;
            filename = "random_scene.png"
        }

        2 => {
            world = Arc::new(two_spheres());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            filename = "two_spheres.png"
        }

        3 => {
            world = Arc::new(two_perlin_spheres());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            filename = "two_perlin_spheres.png"
        }

        4 => {
            world = Arc::new(earth());
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            filename = "earth.png"
        }

        5 => {
            world = Arc::new(simple_light());
            samples_per_pixel = 400;
            background = Vec3::default();
            lookfrom = Vec3::new(26.0, 3.0, 6.0);
            lookat = Vec3::new(0.0, 2.0, 0.0);
            vfov = 20.0;
            filename = "simple_light.png"
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
            filename = "cornell_box.png"
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
            filename = "cornell_box_smoke.png"
        }

        8 | _ => {
            world = Arc::new(the_next_week_final_scene());
            aspect_ratio = 1.0;
            image_width = 800;
            image_height = 800;
            samples_per_pixel = 10000;
            background = Vec3::default();
            lookfrom = Vec3::new(478.0, 278.0, -600.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;
            filename = "the_next_week_final_scene.png"
        }
    }

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

    // Render
    // println!("P3");
    // println!("{} {}", image_width, image_height);
    // println!("255");

    // for j in (0..image_height).rev() {
    //     eprintln!("Scanlines remaining: {}", j);
    //     for i in 0..image_width {
    //         let mut pixel_color = Vec3::new(0.0, 0.0, 0.0);
    //         for _s in 0..samples_per_pixel {
    //             let mut rng = rand::thread_rng();

    //             let u: f64 = (i as f64 + rng.gen::<f64>()) / (image_width - 1) as f64;
    //             let v: f64 = (j as f64 + rng.gen::<f64>()) / (image_height - 1) as f64;
    //             let r: Ray = camera.get_ray(u, v);
    //             pixel_color = pixel_color + ray_color(&r, background, &world, max_depth);
    //         }

    //         let scale: f64 = 1.0 / samples_per_pixel as f64;
    //         println!(
    //             "{} {} {}",
    //             ((256 as f64 * clamp((pixel_color.r() * scale).sqrt(), 0.0, 0.999)) as i32),
    //             ((256 as f64 * clamp((pixel_color.g() * scale).sqrt(), 0.0, 0.999)) as i32),
    //             ((256 as f64 * clamp((pixel_color.b() * scale).sqrt(), 0.0, 0.999)) as i32),
    //         );
    //     }
    // }

    // eprintln!("Done.");

    let n_workers = 8;
    let mut img = image::RgbaImage::new(image_width as u32, image_height as u32);
    let n_jobs = image_width * image_height;
    let pool = ThreadPool::new(n_workers);

    let start_time = time::Instant::now();

    let (tx, rx) = channel();
    for j in (0..image_height) {
        for i in 0..image_width {
            let scale: f64 = 1.0 / samples_per_pixel as f64;
            let tx = tx.clone();
            let camera = camera.clone();
            let world = world.clone();

            pool.execute(move || {
                let mut pixel_color = Vec3::new(0.0, 0.0, 0.0);
                for _s in 0..samples_per_pixel {
                    let mut rng = rand::thread_rng();

                    let u: f64 = (i as f64 + rng.gen::<f64>()) / (image_width - 1) as f64;
                    let v: f64 = 1.0 - (j as f64 + rng.gen::<f64>()) / (image_height - 1) as f64;
                    let r: Ray = camera.get_ray(u, v);
                    pixel_color = pixel_color + ray_color(&r, background, &world, max_depth);
                }

                let pixel = image::Rgba([
                    (256 as f64 * clamp((pixel_color.r() * scale).sqrt(), 0.0, 0.999)) as u8,
                    (256 as f64 * clamp((pixel_color.g() * scale).sqrt(), 0.0, 0.999)) as u8,
                    (256 as f64 * clamp((pixel_color.b() * scale).sqrt(), 0.0, 0.999)) as u8,
                    255,
                ]);

                tx.send(RenderResult { pixel, x: i, y: j }).unwrap();
            })
        }
    }

    for (i, result) in rx.iter().enumerate().take(n_jobs as usize) {
        let result = result as RenderResult;
        img.put_pixel(result.x, result.y, result.pixel)
    }

    eprintln!("{}", (time::Instant::now() - start_time).as_millis() / 1000);

    img.save(format!("{}", filename)).unwrap();
}
