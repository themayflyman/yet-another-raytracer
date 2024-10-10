extern crate rand;

use std::path::Path;
use std::sync::Arc;

use crate::aarect::{XYRect, XZRect, YZRect};
use crate::box_entity::BoxEntity;
use crate::bvh::BVHNode;
use crate::color::RGB;
use crate::hittable::{ConstantMedium, FlipFace, HittableList, RotateY, Translate};
use crate::material::*;
use crate::qbvh::L1QBVH;
use crate::rand::Rng;
use crate::sphere::{MovingSphere, StillSphere};
use crate::texture::{CheckerTexture, ImageTexture, NoiseTexture, NoiseType, SolidColor};
use crate::triangle::{Triangle, TriangleMesh};
use crate::vec3::Vec3;
// use raytracer_codegen::make_static_the_next_week_final_scene;

// make_static_the_next_week_final_scene! {}

pub fn random_scene() -> HittableList {
    let mut world: HittableList = HittableList::new();

    // let checker = CheckerTexture::new(Vec3::new(0.2, 0.3, 0.1), Vec3::new(0.9, 0.9, 0.9));
    let group_material = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
    world.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        group_material,
    )));

    for a in -11..11 {
        for b in -11..11 {
            let mut rng = rand::thread_rng();
            let choose_mat: f32 = rng.gen::<f32>();
            let center: Vec3 = Vec3::new(
                a as f32 + 0.9 * rng.gen::<f32>(),
                0.2,
                b as f32 + 0.9 * rng.gen::<f32>(),
            );

            if (center - Vec3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                if choose_mat < 0.8 {
                    let albedo = SolidColor::new(RGB::new(
                        rng.gen_range(-1.0..1.0),
                        rng.gen_range(-1.0..1.0),
                        rng.gen_range(-1.0..1.0),
                    ));
                    let sphere_material = Lambertian::new(albedo);
                    // let center2 = center + Vec3::new(0.0, rng.gen_range(0.0, 0.5), 0.0);
                    world.add_object(Arc::new(StillSphere::new(
                        center,
                        // center2,
                        // 0.0,
                        // 1.0,
                        0.2,
                        sphere_material,
                    )));
                } else if choose_mat < 0.95 {
                    let albedo = SolidColor::new(RGB::new(
                        rng.gen_range(0.5..1.0),
                        rng.gen_range(0.5..1.0),
                        rng.gen_range(0.5..1.0),
                    ));
                    let fuzz: f32 = rng.gen_range(0.0..0.5);
                    let sphere_material = Metal::new(albedo, fuzz);
                    world.add_object(Arc::new(StillSphere::new(center, 0.2, sphere_material)));
                } else {
                    let sphere_material = SF66;
                    world.add_object(Arc::new(StillSphere::new(center, 0.2, sphere_material)));
                }
            }
        }
    }

    world.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 1.0, 0.0),
        1.0,
        SF66,
    )));

    world.add_object(Arc::new(StillSphere::new(
        Vec3::new(-4.0, 1.0, 0.0),
        1.0,
        Lambertian::new(SolidColor::new(RGB::new(0.4, 0.2, 0.1))),
    )));

    world.add_object(Arc::new(StillSphere::new(
        Vec3::new(4.0, 1.0, 0.0),
        1.0,
        Metal::new(SolidColor::new(RGB::new(0.7, 0.6, 0.5)), 0.0),
    )));

    world
}

pub fn two_spheres() -> HittableList {
    let mut objects = HittableList::new();

    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, -10.0, 0.0),
        10.0,
        Lambertian::new(CheckerTexture::new(
            RGB::new(0.2, 0.3, 0.1),
            RGB::new(0.9, 0.9, 0.9),
        )),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 10.0, 0.0),
        10.0,
        Lambertian::new(CheckerTexture::new(
            RGB::new(0.2, 0.3, 0.1),
            RGB::new(0.9, 0.9, 0.9),
        )),
    )));

    objects
}

pub fn two_perlin_spheres() -> HittableList {
    let mut objects = HittableList::new();

    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 2.0, 0.0),
        2.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
        // Lambertian::new(LinearGradientTexture::new(
        //     vec![ColorStop { color: Vec3::new(1.0, 1.0, 1.0), stop: 0.0 }, ColorStop { color: Vec3::new(0.8, 0.5, 0.3), stop: 0.5 }, ColorStop { color: Vec3::new(0.2, 0.3, 0.1), stop: 1.0 }]
        // )),
    )));

    objects
}

pub fn earth() -> HittableList {
    let mut objects = HittableList::new();

    let earth_texture = ImageTexture::new("input/earthmap.jpg").unwrap();
    let earth_surface = Lambertian::new(earth_texture);
    let global = StillSphere::new(Vec3::new(0.0, 0.0, 0.0), 2.0, earth_surface);
    objects.add_object(Arc::new(global));

    objects
}

pub fn simple_light() -> HittableList {
    let mut objects = HittableList::new();

    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 2.0, 0.0),
        2.0,
        Lambertian::new(NoiseTexture::new(NoiseType::Marble, 4.0)),
    )));

    let difflight = DiffuseLight::new(SolidColor::new(RGB::new(4.0, 4.0, 4.0)));
    objects.add_object(Arc::new(XYRect::new(3.0, 5.0, 1.0, 3.0, -2.0, difflight)));

    objects
}

pub fn cornell_box() -> HittableList {
    let mut objects = HittableList::new();

    let red = Lambertian::new(SolidColor::new(RGB::new(0.65, 0.05, 0.05)));
    let white = Lambertian::new(SolidColor::new(RGB::new(0.73, 0.73, 0.73)));
    let green = Lambertian::new(SolidColor::new(RGB::new(0.12, 0.45, 0.15)));
    let light = DiffuseLight::new(SolidColor::new(RGB::new(15.0, 15.0, 15.0)));

    objects.add_object(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 555.0, green)));
    objects.add_object(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 0.0, red)));
    objects.add_object(Arc::new(FlipFace::new(Arc::new(XZRect::new(
        213.0, 343.0, 227.0, 332.0, 554.0, light,
    )))));
    objects.add_object(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        white.clone(),
    )));
    objects.add_object(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));
    objects.add_object(Arc::new(XYRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));

    let box1 = Arc::new(Translate::new(
        Arc::new(RotateY::new(
            Arc::new(BoxEntity::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(165.0, 330.0, 165.0),
                white,
            )),
            15.0,
        )),
        Vec3::new(265.0, 0.0, 295.0),
    ));
    objects.add_object(box1);

    // let box2 = Arc::new(Translate::new(
    //     Arc::new(RotateY::new(
    //         Arc::new(BoxEntity::new(
    //             Vec3::new(0.0, 0.0, 0.0),
    //             Vec3::new(165.0, 165.0, 165.0),
    //             white,
    //         )),
    //         -18.0,
    //     )),
    //     Vec3::new(130.0, 0.0, 65.0),
    // ));
    // objects.add_object(box2);

    let glass = SF66;
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(190.0, 90.0, 190.0),
        90.0,
        glass,
    )));

    objects
}

pub fn cornell_box_smoke() -> HittableList {
    let mut objects = HittableList::new();

    let red = Lambertian::new(SolidColor::new(RGB::new(0.65, 0.05, 0.05)));
    let white = Lambertian::new(SolidColor::new(RGB::new(0.73, 0.73, 0.73)));
    let green = Lambertian::new(SolidColor::new(RGB::new(0.12, 0.45, 0.15)));
    let light = DiffuseLight::new(SolidColor::new(RGB::new(7.0, 7.0, 7.0)));

    objects.add_object(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 555.0, green)));
    objects.add_object(Arc::new(YZRect::new(0.0, 555.0, 0.0, 555.0, 0.0, red)));
    objects.add_object(Arc::new(FlipFace::new(Arc::new(XZRect::new(
        113.0, 443.0, 127.0, 432.0, 554.0, light,
    )))));
    objects.add_object(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        white.clone(),
    )));
    objects.add_object(Arc::new(XZRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));
    objects.add_object(Arc::new(XYRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));

    let box1 = Arc::new(ConstantMedium::new(
        Translate::new(
            Arc::new(RotateY::new(
                Arc::new(BoxEntity::new(
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(165.0, 330.0, 165.0),
                    white.clone(),
                )),
                15.0,
            )),
            Vec3::new(265.0, 0.0, 295.0),
        ),
        0.01,
        SolidColor::new(RGB::new(0.0, 0.0, 0.0)),
    ));
    let box2 = Arc::new(ConstantMedium::new(
        Translate::new(
            Arc::new(RotateY::new(
                Arc::new(BoxEntity::new(
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(165.0, 165.0, 165.0),
                    white,
                )),
                -18.0,
            )),
            Vec3::new(130.0, 0.0, 65.0),
        ),
        0.01,
        SolidColor::new(RGB::new(1.0, 1.0, 1.0)),
    ));

    objects.add_object(box1);
    objects.add_object(box2);

    objects
}

pub fn the_next_week_final_scene() -> HittableList {
    let mut rng = rand::thread_rng();

    let mut boxes1 = HittableList::new();
    let ground = Lambertian::new(SolidColor::new(RGB::new(0.48, 0.83, 0.53)));

    for i in 0..20 {
        for j in 0..20 {
            let w = 100.0;
            let x0: f32 = -1000.0 + i as f32 * w;
            let z0: f32 = -1000.0 + j as f32 * w;
            let y0: f32 = 0.0;
            let x1: f32 = x0 + w;
            let y1: f32 = rng.gen_range(1.0..101.0);
            let z1: f32 = z0 + w;

            boxes1.add_object(Arc::new(BoxEntity::new(
                Vec3::new(x0, y0, z0),
                Vec3::new(x1, y1, z1),
                ground.clone(),
            )));
        }
    }

    let mut objects = HittableList::new();
    let boxes1_size = boxes1.size();

    objects.add_object(Arc::new(BVHNode::new(
        &mut boxes1.objects,
        0,
        boxes1_size,
        0.0,
        0.0,
    )));

    let light = DiffuseLight::new(SolidColor::new(RGB::new(7.0, 7.0, 7.0)));
    objects.add_object(Arc::new(XZRect::new(
        123.0, 423.0, 147.0, 412.0, 554.0, light,
    )));

    let center0 = Vec3::new(400.0, 400.0, 200.0);
    let center1 = center0 + Vec3::new(30.0, 0.0, 0.0);
    let moving_sphere_material = Lambertian::new(SolidColor::new(RGB::new(0.7, 0.3, 0.1)));
    objects.add_object(Arc::new(MovingSphere::new(
        center0,
        center1,
        0.0,
        1.0,
        50.0,
        moving_sphere_material,
    )));

    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(260.0, 150.0, 45.0),
        50.0,
        SF66,
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 150.0, 145.0),
        50.0,
        Metal::new(SolidColor::new(RGB::new(0.8, 0.8, 0.9)), 1.0),
    )));

    let mut boundary = StillSphere::new(Vec3::new(360.0, 150.0, 145.0), 70.0, SF66);
    let boundary2 = StillSphere::new(Vec3::new(360.0, 150.0, 145.0), 70.0, SF66);
    objects.add_object(Arc::new(boundary));
    objects.add_object(Arc::new(ConstantMedium::new(
        boundary2,
        0.2,
        SolidColor::new(RGB::new(0.2, 0.4, 0.9)),
    )));
    boundary = StillSphere::new(Vec3::new(0.0, 0.0, 0.0), 5000.0, SF66);
    objects.add_object(Arc::new(ConstantMedium::new(
        boundary,
        0.0001,
        SolidColor::new(RGB::new(1.0, 1.0, 1.0)),
    )));

    let emat = Lambertian::new(ImageTexture::new("input/earthmap.jpg").unwrap());
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(400.0, 200.0, 400.0),
        100.0,
        emat,
    )));
    let pertext = NoiseTexture::new(NoiseType::Marble, 0.1);
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(220.0, 280.0, 300.0),
        80.0,
        Lambertian::new(pertext),
    )));

    let mut boxes2 = HittableList::new();
    let boxes2_size = boxes2.size();
    let white = Lambertian::new(SolidColor::new(RGB::new(0.73, 0.73, 0.73)));
    for _ in 0..1000 {
        boxes2.add_object(Arc::new(StillSphere::new(
            Vec3::random(0.0, 165.0),
            10.0,
            white.clone(),
        )));
    }

    objects.add_object(Arc::new(Translate::new(
        Arc::new(RotateY::new(
            Arc::new(BVHNode::new(&mut boxes2.objects, 0, boxes2_size, 0.0, 0.0)),
            15.0,
        )),
        Vec3::new(-100.0, 270.0, 395.0),
    )));

    objects
}

pub fn sycee() -> HittableList {
    let mut objects = HittableList::new();
    let light = DiffuseLight::new(SolidColor::new(RGB::new(5.0, 5.0, 5.0)));
    let ground = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
    let glass = SF66;
    // let metal = Metal::new(SolidColor::new(RGB::new(0.7, 0.6, 0.5)), 0.0);
    // let red = Lambertian::new(SolidColor::new(RGB::new(0.65, 0.05, 0.05)));

    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 6.0, 2.0),
        2.0,
        light,
    )));
    objects.add_object(Arc::new(TriangleMesh::from_obj(
        Path::new("input/sycee.obj"),
        glass,
    )));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-20.0, 0.0, -30.0),
            Vec3::new(20.0, 0.0, -30.0),
            Vec3::new(20.0, 0.0, 30.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        material: ground.clone(),
    }));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-20.0, 0.0, -30.0),
            Vec3::new(-20.0, 0.0, 30.0),
            Vec3::new(20.0, 0.0, 30.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        material: ground,
    }));

    return objects;
}

pub fn teapot() -> HittableList {
    let mut objects = HittableList::new();
    let light = DiffuseLight::new(SolidColor::new(RGB::new(5.0, 5.0, 5.0)));
    let ground = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
    let glass = SF66;
    // let silver = Metal::new(SolidColor::new(RGB::new(0.7976, 0.7976, 0.7976)), 0.0);

    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(30.0, 40.0, -30.0),
        20.0,
        light.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(-20.0, 10.0, 50.0),
        10.0,
        light,
    )));
    objects.add_object(Arc::new(TriangleMesh::from_obj(
        Path::new("input/teapot.obj"),
        glass,
    )));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-80.0, 0.0, -120.0),
            Vec3::new(80.0, 0.0, -120.0),
            Vec3::new(80.0, 0.0, 120.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        material: ground.clone(),
    }));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-80.0, 0.0, -120.0),
            Vec3::new(-80.0, 0.0, 120.0),
            Vec3::new(80.0, 0.0, 120.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        material: ground,
    }));

    return objects;
}

pub fn bunny() -> HittableList {
    let mut objects = HittableList::new();
    let glass = SF66;
    let ground = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
    let light = DiffuseLight::new(SolidColor::new(RGB::new(5.0, 5.0, 5.0)));
    objects.add_object(Arc::new(TriangleMesh::from_obj(
        Path::new("input/bunny.obj"),
        glass,
    )));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-20.0, 0.0, -30.0),
            Vec3::new(20.0, 0.0, -30.0),
            Vec3::new(20.0, 0.0, 30.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        material: ground.clone(),
    }));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-20.0, 0.0, -30.0),
            Vec3::new(-20.0, 0.0, 30.0),
            Vec3::new(20.0, 0.0, 30.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        material: ground,
    }));
    objects.add_object(Arc::new(StillSphere::new(
Vec3::new(0.0, 6.0, -2.0),
        2.0,
        light.clone(),
    )));

    return objects;
}

pub fn david() -> HittableList {
    let mut objects = HittableList::new();
    let light = DiffuseLight::new(SolidColor::new(RGB::new(5.0, 5.0, 5.0)));
    let white = Lambertian::new(SolidColor::new(RGB::new(1.0, 1.0, 1.0)));
    objects.add_object(Arc::new(TriangleMesh::from_obj(
        Path::new("input/david.obj"),
        white,
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(1200.0, 1300.0, 800.0),
        700.0,
        light.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(-1200.0, 1300.0, 800.0),
        700.0,
        light.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(1200.0, 1300.0, -800.0),
        700.0,
        light.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(1200.0, -1300.0, -800.0),
        700.0,
        light.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(1200.0, 1300.0, -800.0),
        700.0,
        light,
    )));

    return objects;
}

pub fn three_spheres() -> HittableList {
    let glass = SF66;
    let ground = Lambertian::new(SolidColor::new(RGB::new(0.5, 0.5, 0.5)));
    let light = DiffuseLight::new(SolidColor::new(RGB::new(5.0, 5.0, 5.0)));
    // let earth = Lambertian::new(ImageTexture::new("input/earthmap.jpg").unwrap());
    // let metal = Metal::new(SolidColor::new(RGB::new(0.7, 0.6, 0.5)), 0.0);

    let mut objects = HittableList::new();
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-20.0, 0.0, -30.0),
            Vec3::new(20.0, 0.0, -30.0),
            Vec3::new(20.0, 0.0, 30.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        material: ground.clone(),
    }));
    objects.add_object(Arc::new(Triangle {
        vertices: [
            Vec3::new(-20.0, 0.0, -30.0),
            Vec3::new(-20.0, 0.0, 30.0),
            Vec3::new(20.0, 0.0, 30.0),
        ],
        normals: [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        uv: [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        material: ground,
    }));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 1.0, 0.0),
        1.0,
        glass,
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 1.3, 0.0),
        -0.7,
        glass.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 0.65, 0.0),
        -0.35,
        glass.clone(),
    )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 0.325, 0.0),
        -0.125,
        glass.clone(),
    )));
    // objects.add_object(Arc::new(StillSphere::new(
    //     Vec3::new(-3.0, 1.0, 0.0),
    //     1.0,
    //     earth,
    // )));
    // objects.add_object(Arc::new(StillSphere::new(
    //     Vec3::new(3.0, 1.0, 0.0),
    //     1.0,
    //     metal,
    // )));
    objects.add_object(Arc::new(StillSphere::new(
        Vec3::new(0.0, 6.0, 2.0),
        2.0,
        light,
    )));

    return objects;
}
