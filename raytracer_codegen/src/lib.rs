extern crate proc_macro;
use proc_macro2::TokenStream;
use quote::quote;

use rand::Rng;

#[derive(Copy, Clone, Debug)]
struct Vec3(f64, f64, f64);

#[derive(Copy, Clone, Debug)]
struct StillSphere {
    pub pos: Vec3,
    pub size: f64,
    pub color: Vec3,
}

#[derive(Copy, Clone, Debug)]
struct BoxEntity {
    pub p0: Vec3,
    pub p1: Vec3,
}

fn generate_boxes() -> Vec<BoxEntity> {
    let mut rng = rand::thread_rng();
    let mut boxes = vec![];
    for i in 0..20 {
        for j in 0..20 {
            let w = 100.0;
            let x0: f64 = -1000.0 + i as f64 * w;
            let z0: f64 = -1000.0 + j as f64 * w;
            let y0: f64 = 0.0;
            let x1: f64 = x0 + w;
            let y1: f64 = rng.gen_range(1.0, 101.0);
            let z1: f64 = z0 + w;

            boxes.push(BoxEntity {
                p0: Vec3(x0, y0, z0),
                p1: Vec3(x1, y1, z1),
            });
        }
    }
    boxes
}

fn generate_spheres() -> Vec<StillSphere> {
    let mut rng = rand::thread_rng();
    let mut spheres = vec![];
    for _ in 0..1000 {
        let pos = Vec3 (
            rng.gen_range(0.0, 165.0),
            rng.gen_range(0.0, 165.0),
            rng.gen_range(0.0, 165.0),
        );
        let size = 10.0;
        spheres.push(StillSphere {
            pos,
            size,
            color: Vec3 ( 0.73, 0.73, 0.73 )
        });
    }
    spheres
}

fn sphere_code(sphere: &StillSphere) -> TokenStream {
    let Vec3(x, y, z) = sphere.pos;
    let Vec3(r, g, b) = sphere.color;
    let size = sphere.size * 0.9;

    quote! {
        Arc::new(StillSphere::new(
            Vec3::new(#x, #y, #z),
            #size,
            Lambertian::new(SolidColor::new_from_value(#r, #g, #b))
        ))
    }
}

fn box_code(b: &BoxEntity) -> TokenStream {
    let Vec3(x0, y0, z0) = b.p0;
    let Vec3(x1, y1, z1) = b.p1;

    quote! {
        Arc::new(BoxEntity::new(
            Vec3::new(#x0, #y0, #z0),
            Vec3::new(#x1, #y1, #z1),
            Lambertian::new(SolidColor::new(Vec3::new(0.48, 0.83, 0.53))),
        ))
    }
}

fn construct_box_bvh(mut boxes: Vec<BoxEntity>) -> TokenStream {
    match boxes.len() {
        0 => unreachable!(),
        1 => {
            let box_left = box_code(boxes.first().unwrap());
            let box_right = box_code(boxes.first().unwrap());
            quote! { Arc::new(BVHNodeStatic::construct(#box_left, #box_right)) }
        }
        _ => {
            let axis = rand::thread_rng().gen_range(0, 3);
            boxes.sort_by(|a, b| match axis {
                0 => a.p0.0.partial_cmp(&b.p0.0).unwrap(),
                1 => a.p0.1.partial_cmp(&b.p0.1).unwrap(),
                2 => a.p0.2.partial_cmp(&b.p0.2).unwrap(),
                _ => unreachable!(),
            });
            let mut a = boxes;
            let b = a.split_off(a.len() / 2);
            let left = construct_box_bvh(a);
            let right = construct_box_bvh(b);
            quote! { Arc::new(BVHNodeStatic::construct(#left, #right)) }
        }
    }
}

fn construct_sphere_bvh(mut spheres: Vec<StillSphere>) -> TokenStream {
    match spheres.len() {
        0 => unreachable!(),
        1 => {
            let sphere_left = sphere_code(spheres.first().unwrap());
            let sphere_right = sphere_code(spheres.first().unwrap());
            quote! { Arc::new(BVHNodeStatic::construct(#sphere_left, #sphere_right)) }
        }
        _ => {
            let axis = rand::thread_rng().gen_range(0, 3);
            spheres.sort_by(|a, b| match axis {
                0 => a.pos.0.partial_cmp(&b.pos.0).unwrap(),
                1 => a.pos.1.partial_cmp(&b.pos.1).unwrap(),
                2 => a.pos.2.partial_cmp(&b.pos.2).unwrap(),
                _ => unreachable!(),
            });
            let mut a = spheres;
            let b = a.split_off(a.len() / 2);
            let left = construct_sphere_bvh(a);
            let right = construct_sphere_bvh(b);
            quote! { Arc::new(BVHNodeStatic::construct(#left, #right)) }
        }
    }
}

#[proc_macro]
pub fn make_static_the_next_week_final_scene(_item: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let boxes = generate_boxes();
    let boxes_bvh: TokenStream = construct_box_bvh(boxes);
    let spheres = generate_spheres();
    let spheres_bvh: TokenStream = construct_sphere_bvh(spheres);
    proc_macro::TokenStream::from(quote! {
        pub fn static_the_next_week_final_scene() -> HittableList {
            let mut objects = HittableList::new();

            let light = DiffuseLight::new(SolidColor::new(Vec3::new(7.0, 7.0, 7.0)));
            objects.add_object(Arc::new(FlipFace::new(Arc::new(XZRect::new(
                123.0, 423.0, 147.0, 412.0, 554.0, light,
            )))));

            objects.add_object(#boxes_bvh);

            let center0 = Vec3::new(400.0, 400.0, 200.0);
            let center1 = center0 + Vec3::new(30.0, 0.0, 0.0);
            let moving_sphere_material = Lambertian::new(SolidColor::new(Vec3::new(0.7, 0.3, 0.1)));
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
                Dielectric::new(1.5),
            )));
            objects.add_object(Arc::new(StillSphere::new(
                Vec3::new(0.0, 150.0, 145.0),
                50.0,
                Metal::new(Vec3::new(0.8, 0.8, 0.9), 1.0),
            )));

            let mut boundary = StillSphere::new(Vec3::new(360.0, 150.0, 145.0), 70.0, Dielectric::new(1.5));
            let boundary2 = StillSphere::new(Vec3::new(360.0, 150.0, 145.0), 70.0, Dielectric::new(1.5));
            objects.add_object(Arc::new(boundary));
            objects.add_object(Arc::new(ConstantMedium::new(
                boundary2,
                0.2,
                SolidColor::new(Vec3::new(0.2, 0.4, 0.9)),
            )));
            boundary = StillSphere::new(Vec3::new(0.0, 0.0, 0.0), 5000.0, Dielectric::new(1.5));
            objects.add_object(Arc::new(ConstantMedium::new(
                boundary,
                0.0001,
                SolidColor::new(Vec3::new(1.0, 1.0, 1.0)),
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

            objects.add_object(Arc::new(Translate::new(
                Arc::new(RotateY::new(
                    #spheres_bvh,
                    15.0,
                )),
                Vec3::new(-100.0, 270.0, 395.0),
            )));

            objects
        }
    })
}
