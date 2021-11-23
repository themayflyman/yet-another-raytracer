extern crate rand;

use camera::Camera;
use hittable::{Hittable, HittableList};
use material::{Dielectric, Lambertian, Material, Metal};
use rand::Rng;
use ray::Ray;
use sphere::{MovingSphere, StillSphere};
use texture::{CheckerTexture, NoiseTexture, SolidColor};
use vec3::Vec3;

use self::aarect::{XYRect, XZRect, YZRect};
use self::material::DiffuseLight;
use self::texture::{ImageTexture, NoiseType};

mod aabb;
mod aarect;
mod bvh;
mod camera;
mod hittable;
mod material;
mod ray;
mod sphere;
mod texture;
mod vec3;

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
        let emitted = rec.material.emitted(rec.u, rec.v, rec.p);

        if let Some(scattered_ray) = scattered.ray {
            return emitted + scattered.color * ray_color(&scattered_ray, background, world, depth - 1);
        } else {
            return emitted;
        }
    }

    return background;
}

fn random_scene() -> HittableList {
    let mut world: HittableList = HittableList::new();

    let checker = CheckerTexture::new(Vec3::new(0.2, 0.3, 0.1), Vec3::new(0.9, 0.9, 0.9));
    let group_material = Box::new(Lambertian::new(Box::new(checker)));
    world.add_sphere(Box::new(StillSphere::new(
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
                let sphere_material: Box<dyn Material>;

                if choose_mat < 0.8 {
                    let albedo = Box::new(SolidColor::new(Vec3::new(
                        rng.gen_range(-1.0, 1.0),
                        rng.gen_range(-1.0, 1.0),
                        rng.gen_range(-1.0, 1.0),
                    )));
                    sphere_material = Box::new(Lambertian::new(albedo));
                    let center2 = center + Vec3::new(0.0, rng.gen_range(0.0, 0.5), 0.0);
                    world.add_sphere(Box::new(MovingSphere::new(
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
                    sphere_material = Box::new(Metal::new(albedo, fuzz));
                    world.add_sphere(Box::new(StillSphere::new(center, 0.2, sphere_material)));
                } else {
                    sphere_material = Box::new(Dielectric::new(1.5));
                    world.add_sphere(Box::new(StillSphere::new(center, 0.2, sphere_material)));
                }
            }
        }
    }

    let material1 = Box::new(Dielectric::new(1.5));
    world.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, 1.0, 0.0),
        1.0,
        material1,
    )));

    let material2 = Box::new(Lambertian::new(Box::new(SolidColor::new(Vec3::new(
        0.4, 0.2, 0.1,
    )))));
    world.add_sphere(Box::new(StillSphere::new(
        Vec3::new(-4.0, 1.0, 0.0),
        1.0,
        material2,
    )));

    let material3 = Box::new(Metal::new(Vec3::new(0.7, 0.6, 0.5), 0.0));
    world.add_sphere(Box::new(StillSphere::new(
        Vec3::new(4.0, 1.0, 0.0),
        1.0,
        material3,
    )));

    return world;
}

fn two_spheres() -> HittableList {
    let mut objects = HittableList::new();

    let checker = CheckerTexture::new(Vec3::new(0.2, 0.3, 0.1), Vec3::new(0.9, 0.9, 0.9));
    let checker2 = checker.clone();

    objects.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, -10.0, 0.0),
        10.0,
        Box::new(Lambertian::new(Box::new(checker))),
    )));
    objects.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, 10.0, 0.0),
        10.0,
        Box::new(Lambertian::new(Box::new(checker2))),
    )));

    objects
}

fn two_perlin_spheres() -> HittableList {
    let mut objects = HittableList::new();

    let pertext = NoiseTexture::new(NoiseType::Marble, 4.0);
    let pertext2 = pertext.clone();

    objects.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Box::new(Lambertian::new(Box::new(pertext))),
    )));
    objects.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, 2.0, 0.0),
        2.0,
        Box::new(Lambertian::new(Box::new(pertext2))),
    )));

    objects
}

fn earth() -> HittableList {
    let mut objects = HittableList::new();

    let earth_texture = ImageTexture::new("earthmap.jpg").unwrap();
    let earth_surface = Lambertian::new(Box::new(earth_texture));
    let global = StillSphere::new(Vec3::new(0.0, 0.0, 0.0), 2.0, Box::new(earth_surface));
    objects.add_sphere(Box::new(global));

    objects
}

fn simple_light() -> HittableList {
    let mut objects = HittableList::new();

    let pertext = NoiseTexture::new(NoiseType::Marble, 4.0);
    let pertext2 = pertext.clone();
    objects.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Box::new(Lambertian::new(Box::new(pertext))),
    )));
    objects.add_sphere(Box::new(StillSphere::new(
        Vec3::new(0.0, 2.0, 0.0),
        2.0,
        Box::new(Lambertian::new(Box::new(pertext2))),
    )));

    let difflight = DiffuseLight::new(Box::new(SolidColor::new(Vec3::new(4.0, 4.0, 4.0))));
    objects.add_sphere(Box::new(XYRect::new(
        3.0,
        5.0,
        1.0,
        3.0,
        -2.0,
        Box::new(difflight),
    )));

    objects
}

fn cornell_box() -> HittableList {
    let mut objects = HittableList::new();

    let red = Box::new(Lambertian::new(Box::new(SolidColor::new(Vec3::new(
        0.65, 0.05, 0.05,
    )))));
    let white = Box::new(Lambertian::new(Box::new(SolidColor::new(Vec3::new(
        0.73, 0.73, 0.73,
    )))));
    let green = Box::new(Lambertian::new(Box::new(SolidColor::new(Vec3::new(
        0.12, 0.45, 0.15,
    )))));
    let light = Box::new(DiffuseLight::new(Box::new(SolidColor::new(Vec3::new(
        15.0, 15.0, 15.0,
    )))));

    objects.add_sphere(Box::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 555.0, green)));
    objects.add_sphere(Box::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 0.0, red)));
    objects.add_sphere(Box::new(XZRect::new(
        213.0, 343.0, 227.0, 332.0, 554.0, light,
    )));
    objects.add_sphere(Box::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        white.clone(),
    )));
    objects.add_sphere(Box::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));
    objects.add_sphere(Box::new(XYRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));

    objects
}

fn main() {
    // Image
    let mut aspect_ratio: f64 = 16.0 / 9.0;
    let mut image_width: usize = 400;
    let mut image_height: usize = (image_width as f64 / aspect_ratio) as usize;
    let mut max_depth: usize = 50;

    // World
    let world: HittableList;
    let lookfrom: Vec3;
    let lookat: Vec3;
    let vfov: f64;
    let mut aperture = 0.0;
    let background: Vec3;
    let mut samples_per_pixel: usize = 100;

    let scene = 6;

    match scene {
        1 => {
            world = random_scene();
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
            aperture = 0.1;
        }

        2 => {
            world = two_spheres();
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
        }

        3 => {
            world = two_perlin_spheres();
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
        }

        4 => {
            world = earth();
            background = Vec3::new(0.7, 0.8, 1.0);
            lookfrom = Vec3::new(13.0, 2.0, 3.0);
            lookat = Vec3::new(0.0, 0.0, 0.0);
            vfov = 20.0;
        }

        5 => {
            world = simple_light();
            samples_per_pixel = 400;
            background = Vec3::default();
            lookfrom = Vec3::new(26.0, 3.0, 6.0);
            lookat = Vec3::new(0.0, 2.0, 0.0);
            vfov = 20.0;
        }

        6 | _ => {
            world = cornell_box();
            aspect_ratio = 1.0;
            image_width = 600;
            image_height = 600;
            samples_per_pixel = 200;
            background = Vec3::default();
            lookfrom = Vec3::new(278.0, 278.0, -800.0);
            lookat = Vec3::new(278.0, 278.0, 0.0);
            vfov = 40.0;
        }
    }

    let vup: Vec3 = Vec3::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;

    // Camera
    let camera: Camera = Camera::new(
        lookfrom,
        lookat,
        vup,
        vfov,
        aspect_ratio,
        aperture,
        dist_to_focus,
        0.0,
        1.0,
    );

    // Render
    println!("P3");
    println!("{} {}", image_width, image_height);
    println!("255");

    for j in (0..image_height).rev() {
        eprintln!("Scanlines remaining: {}", j);
        for i in 0..image_width {
            let mut pixel_color = Vec3::new(0.0, 0.0, 0.0);
            for _s in 0..samples_per_pixel {
                let mut rng = rand::thread_rng();

                let u: f64 = (i as f64 + rng.gen::<f64>()) / (image_width - 1) as f64;
                let v: f64 = (j as f64 + rng.gen::<f64>()) / (image_height - 1) as f64;
                let r: Ray = camera.get_ray(u, v);
                pixel_color = pixel_color + ray_color(&r, background, &world, max_depth);
            }

            let scale: f64 = 1.0 / samples_per_pixel as f64;
            println!(
                "{} {} {}",
                ((256 as f64 * clamp((pixel_color.r() * scale).sqrt(), 0.0, 0.999)) as i32),
                ((256 as f64 * clamp((pixel_color.g() * scale).sqrt(), 0.0, 0.999)) as i32),
                ((256 as f64 * clamp((pixel_color.b() * scale).sqrt(), 0.0, 0.999)) as i32),
            );
        }
    }

    eprintln!("Done.");
}
