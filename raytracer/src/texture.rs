use std::convert::TryInto;
use std::error::Error;
// use std::f32::consts::PI;
use std::path::Path;

// use crate::clamp;
use crate::color::{HasReflectance, RGB};
use crate::hittable::HitRecord;
use crate::ray::Ray;
use crate::vec3::Vec3;
use rand::Rng;

pub trait Texture: Clone + Send + Sync {
    // fn value(&self, u: f32, v: f32, p: Vec3) -> RGB;
    fn value(&self, ray_in: &Ray, hit_record: &HitRecord) -> f32;
}

#[derive(Clone)]
pub struct SolidColor<T: HasReflectance> {
    color_value: T,
}

impl<T: HasReflectance> SolidColor<T> {
    pub fn new(color_value: T) -> Self {
        Self { color_value }
    }

    // #[allow(dead_code)]
    // pub fn new_from_value(red: f32, green: f32, blue: f32) -> Self {
    //     Self {
    //         color_value: RGB::new(red, green, blue),
    //     }
    // }
}

impl<T: HasReflectance> Texture for SolidColor<T> {
    fn value(&self, ray_in: &Ray, _hit_record: &HitRecord) -> f32 {
        return self.color_value.reflect(ray_in.wavelength);
    }
}

#[derive(Clone)]
pub struct CheckerTexture<T: HasReflectance> {
    odd: SolidColor<T>,
    even: SolidColor<T>,
}

impl<T: HasReflectance> CheckerTexture<T> {
    pub fn new(color1: T, color2: T) -> Self {
        Self {
            odd: SolidColor::new(color1),
            even: SolidColor::new(color2),
        }
    }
}

impl<T: HasReflectance> Texture for CheckerTexture<T> {
    fn value(&self, ray_in: &Ray, hit_record: &HitRecord) -> f32 {
        let sines = f32::sin(10.0 * hit_record.p.x())
            * f32::sin(10.0 * hit_record.p.y())
            * f32::sin(10.0 * hit_record.p.z());
        if sines < 0.0 {
            self.odd.value(ray_in, hit_record)
        } else {
            self.even.value(ray_in, hit_record)
        }
    }
}

#[derive(Clone)]
pub enum NoiseType {
    // Unfiltered noise
    Square,
    // Trillinear interpolation
    Trillinear,
    // Smooth interpolation with random unit vectors
    Smooth,
    // Marble pattern with adjustable phase
    Marble,
    // Turbulent pattern that resembles a net
    Net,
}
//
#[derive(Clone)]
pub struct Perlin {
    ranfloat: Vec<f32>,
    ranvec: Vec<Vec3>,
    perm_x: Vec<i32>,
    perm_y: Vec<i32>,
    perm_z: Vec<i32>,
    noise_type: NoiseType,
}

impl Perlin {
    const POINT_COUNT: usize = 256;

    pub fn new(noise_type: NoiseType) -> Self {
        let mut rng = rand::thread_rng();

        Self {
            ranfloat: (0..Self::POINT_COUNT)
                .map(|_| rng.gen_range::<f32>(0.0, 1.0))
                .collect(),
            ranvec: (0..Self::POINT_COUNT)
                .map(|_| Vec3::random(-1.0, 1.0))
                .collect(),
            perm_x: Self::perlin_generate_perm(),
            perm_y: Self::perlin_generate_perm(),
            perm_z: Self::perlin_generate_perm(),
            noise_type,
        }
    }

    pub fn noise(&self, p: Vec3) -> f32 {
        match self.noise_type {
            NoiseType::Square => {
                let i: usize = ((4.0 * p.x()) as i32 & 255).try_into().unwrap();
                let j: usize = ((4.0 * p.y()) as i32 & 255).try_into().unwrap();
                let k: usize = ((4.0 * p.z()) as i32 & 255).try_into().unwrap();

                self.ranfloat[(self.perm_x[i] ^ self.perm_y[j] ^ self.perm_z[k]) as usize]
            }

            NoiseType::Trillinear => {
                let mut u = p.x() - p.x().floor();
                let mut v = p.y() - p.y().floor();
                let mut w = p.z() - p.z().floor();
                u = u * u * (3.0 - 2.0 * u);
                v = v * v * (3.0 - 2.0 * v);
                w = w * w * (3.0 - 2.0 * w);

                let i = p.x().floor() as i32;
                let j = p.y().floor() as i32;
                let k = p.z().floor() as i32;

                let mut c = [[[0.0f32; 2]; 2]; 2];
                for (di, i_item) in c.iter_mut().enumerate() {
                    for (dj, ij_item) in i_item.iter_mut().enumerate() {
                        for (dk, ijk_item) in ij_item.iter_mut().enumerate() {
                            *ijk_item = self.ranfloat[(self.perm_x
                                [((i + di as i32) & 255) as usize]
                                ^ self.perm_y[((j + dj as i32) & 255) as usize]
                                ^ self.perm_z[((k + dk as i32) & 255) as usize])
                                as usize]
                        }
                    }
                }

                Self::trilinear_interp(c, u, v, w)
            }

            NoiseType::Smooth | NoiseType::Marble | NoiseType::Net => {
                let u = p.x() - p.x().floor();
                let v = p.y() - p.y().floor();
                let w = p.z() - p.z().floor();

                let i = p.x().floor() as i32;
                let j = p.y().floor() as i32;
                let k = p.z().floor() as i32;

                let mut c = [[[Vec3::default(); 2]; 2]; 2];
                for (di, i_item) in c.iter_mut().enumerate() {
                    for (dj, ij_item) in i_item.iter_mut().enumerate() {
                        for (dk, ijk_item) in ij_item.iter_mut().enumerate() {
                            *ijk_item = self.ranvec[(self.perm_x[((i + di as i32) & 255) as usize]
                                ^ self.perm_y[((j + dj as i32) & 255) as usize]
                                ^ self.perm_z[((k + dk as i32) & 255) as usize])
                                as usize]
                        }
                    }
                }

                Self::perlin_interp(c, u, v, w)
            }
        }
    }

    fn perlin_generate_perm() -> Vec<i32> {
        let mut p: Vec<i32> = (0..Self::POINT_COUNT).map(|i| i as i32).collect();

        Self::permute(&mut p);

        p
    }

    fn permute(p: &mut Vec<i32>) {
        let mut rnd = rand::thread_rng();

        for i in (1..p.len()).rev() {
            let target = rnd.gen_range::<usize>(0, i);
            p.swap(i, target);
        }
    }

    fn trilinear_interp(c: [[[f32; 2]; 2]; 2], u: f32, v: f32, w: f32) -> f32 {
        let mut accum = 0.0;
        for (i, i_item) in c.iter().enumerate() {
            for (j, ij_item) in i_item.iter().enumerate() {
                for (k, ijk_item) in ij_item.iter().enumerate() {
                    accum += (i as f32 * u + (1 - i) as f32 * (1.0 - u))
                        * (j as f32 * v + (1 - j) as f32 * (1.0 - v))
                        * (k as f32 * w + (1 - k) as f32 * (1.0 - w))
                        * ijk_item;
                }
            }
        }

        accum
    }

    fn perlin_interp(c: [[[Vec3; 2]; 2]; 2], u: f32, v: f32, w: f32) -> f32 {
        let uu = u * u * (3.0 - 2.0 * u);
        let vv = v * v * (3.0 - 2.0 * v);
        let ww = w * w * (3.0 - 2.0 * w);
        let mut accum = 0.0;

        for (i, i_item) in c.iter().enumerate() {
            for (j, ij_item) in i_item.iter().enumerate() {
                for (k, ijk_item) in ij_item.iter().enumerate() {
                    let weight_v = Vec3::new(u - i as f32, v - j as f32, w - k as f32);
                    accum += (i as f32 * uu + (1.0 - i as f32) * (1.0 - uu))
                        * (j as f32 * vv + (1.0 - j as f32) * (1.0 - vv))
                        * (k as f32 * ww + (1.0 - k as f32) * (1.0 - ww))
                        * weight_v.dot(ijk_item);
                }
            }
        }

        accum
    }

    pub fn turb(&self, p: Vec3, depth: usize) -> f32 {
        let mut accum = 0.0;
        let mut temp_p = p;
        let mut weight = 1.0;

        for _ in 0..depth {
            accum += weight * self.noise(temp_p);
            weight *= 0.5;
            temp_p = temp_p * 2.0;
        }

        accum.abs()
    }
}

#[derive(Clone)]
pub struct NoiseTexture {
    noise: Perlin,
    scale: f32,
}

impl NoiseTexture {
    pub fn new(noise_type: NoiseType, scale: f32) -> Self {
        Self {
            noise: Perlin::new(noise_type),
            scale,
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, ray_in: &Ray, hit_record: &HitRecord) -> f32 {
        match self.noise.noise_type {
            NoiseType::Net => {
                RGB::new(1.0, 1.0, 1.0).reflect(ray_in.wavelength)
                    * self.noise.turb(hit_record.p * self.scale, 7)
            }

            NoiseType::Marble => {
                RGB::new(1.0, 1.0, 1.0).reflect(ray_in.wavelength)
                    * 0.5
                    * (1.0
                        + f32::sin(
                            self.scale * hit_record.p.z() + 10.0 * self.noise.turb(hit_record.p, 7),
                        ))
            }

            _ => {
                RGB::new(1.0, 1.0, 1.0).reflect(ray_in.wavelength)
                    * 0.5
                    * (1.0 + self.noise.noise(hit_record.p * self.scale))
            }
        }
    }
}

#[derive(Clone)]
pub struct ImageTexture {
    // Buffers of pixels
    data: Vec<u8>,
    width: u32,
    height: u32,
    bytes_per_pixel: u32,
    bytes_per_scanline: u32,
}

impl ImageTexture {
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<Self, Box<dyn Error>> {
        let bytes_per_pixel = 3;
        let img = image::open(filename)?.to_rgb8();
        let (width, height) = img.dimensions();
        let data = img.into_raw();

        Ok(Self {
            data,
            width,
            height,
            bytes_per_pixel,
            bytes_per_scanline: bytes_per_pixel * width,
        })
    }
}

impl Texture for ImageTexture {
    fn value(&self, ray_in: &Ray, hit_record: &HitRecord) -> f32 {
        // If we have no texture data, then return solid cyan as a debugging aid
        if self.data.is_empty() {
            return 1.0;
        }

        // Clamp input texture coordinates to [0,1] x [1,0]
        let uu = hit_record.u.clamp(0.0, 1.0);
        let vv = 1.0 - hit_record.v.clamp(0.0, 1.0); // Filp V to image coordinates

        let mut i = uu as u32 * self.width;
        let mut j = vv as u32 * self.height;

        // Clamp intefer mapping, since actual coordinates should be less than 1.0
        if i >= self.width {
            i = self.width - 1;
        }
        if j >= self.height {
            j = self.height - 1;
        }

        let color_scale = 1.0 / 255.0;
        let pixel = (j * self.bytes_per_scanline + i * self.bytes_per_pixel) as usize;

        let rgb = RGB::new(
            color_scale * self.data[pixel] as f32,
            color_scale * self.data[pixel + 1] as f32,
            color_scale * self.data[pixel + 2] as f32,
        );
        return rgb.reflect(ray_in.wavelength);
    }
}

// #[derive(Clone, Copy)]
// pub struct ColorStop {
//     pub color: RGB,
//     pub stop: f32,
// }
//
// #[derive(Clone)]
// pub struct LinearGradientTexture {
//     color_stops: Vec<ColorStop>,
// }
//
// impl LinearGradientTexture {
//     pub fn new(color_stops: Vec<ColorStop>) -> Self {
//         Self {
//             color_stops,
//         }
//     }
// }
//
// impl Texture for LinearGradientTexture {
//     fn value(&self, u: f32, v: f32, p: Vec3) -> RGB {
//         let percent = 1.0 - (1.0 - (v * PI).cos()) / 2.0;
//         let end_color_stop_idx: usize = self.color_stops.iter().position(|&x| x.stop > percent).unwrap();
//         let start_color_stop_idx = end_color_stop_idx - 1;
//         let end_color_stop = self.color_stops[end_color_stop_idx];
//         let start_color_stop = self.color_stops[start_color_stop_idx];
//
//         return start_color_stop.color + (percent - start_color_stop.stop) / (end_color_stop.stop - start_color_stop.stop) * (end_color_stop.color - start_color_stop.color);
//     }
// }
