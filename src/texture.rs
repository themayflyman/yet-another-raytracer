use std::convert::TryInto;
use std::error::Error;
use std::path::Path;

use crate::clamp;
use crate::vec3::Vec3;
use rand::Rng;

pub trait Texture: Clone + Send + Sync {
    fn value(&self, u: f64, v: f64, p: Vec3) -> Vec3;
}

#[derive(Clone)]
pub struct SolidColor {
    color_value: Vec3,
}

impl SolidColor {
    pub fn new(color_value: Vec3) -> Self {
        Self { color_value }
    }

    #[allow(dead_code)]
    pub fn new_from_value(red: f64, green: f64, blue: f64) -> Self {
        Self {
            color_value: Vec3::new(red, green, blue),
        }
    }
}

impl Texture for SolidColor {
    fn value(&self, _u: f64, _v: f64, _p: Vec3) -> Vec3 {
        self.color_value
    }
}

#[derive(Clone)]
pub struct CheckerTexture {
    odd: SolidColor,
    even: SolidColor,
}

impl CheckerTexture {
    pub fn new(color1: Vec3, color2: Vec3) -> Self {
        Self {
            odd: SolidColor::new(color1),
            even: SolidColor::new(color2),
        }
    }
}

impl Texture for CheckerTexture {
    fn value(&self, u: f64, v: f64, p: Vec3) -> Vec3 {
        let sines = f64::sin(10.0 * p.x()) * f64::sin(10.0 * p.y()) * f64::sin(10.0 * p.z());
        if sines < 0.0 {
            self.odd.value(u, v, p)
        } else {
            self.even.value(u, v, p)
        }
    }
}

#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
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

#[derive(Clone)]
pub struct Perlin {
    ranfloat: Vec<f64>,
    ranvec: Vec<Vec3>,
    perm_x: Vec<i64>,
    perm_y: Vec<i64>,
    perm_z: Vec<i64>,
    noise_type: NoiseType,
}

impl Perlin {
    const POINT_COUNT: usize = 256;

    pub fn new(noise_type: NoiseType) -> Self {
        let mut rng = rand::thread_rng();

        Self {
            ranfloat: (0..Self::POINT_COUNT)
                .map(|_| rng.gen_range::<f64>(0.0, 1.0))
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

    pub fn noise(&self, p: Vec3) -> f64 {
        match self.noise_type {
            NoiseType::Square => {
                let i: usize = ((4.0 * p.x()) as i64 & 255).try_into().unwrap();
                let j: usize = ((4.0 * p.y()) as i64 & 255).try_into().unwrap();
                let k: usize = ((4.0 * p.z()) as i64 & 255).try_into().unwrap();

                self.ranfloat[(self.perm_x[i] ^ self.perm_y[j] ^ self.perm_z[k]) as usize]
            }

            NoiseType::Trillinear => {
                let mut u = p.x() - p.x().floor();
                let mut v = p.y() - p.y().floor();
                let mut w = p.z() - p.z().floor();
                u = u * u * (3.0 - 2.0 * u);
                v = v * v * (3.0 - 2.0 * v);
                w = w * w * (3.0 - 2.0 * w);

                let i = p.x().floor() as i64;
                let j = p.y().floor() as i64;
                let k = p.z().floor() as i64;

                let mut c = [[[0.0f64; 2]; 2]; 2];
                for (di, i_item) in c.iter_mut().enumerate() {
                    for (dj, ij_item) in i_item.iter_mut().enumerate() {
                        for (dk, ijk_item) in ij_item.iter_mut().enumerate() {
                            *ijk_item = self.ranfloat[(self.perm_x
                                [((i + di as i64) & 255) as usize]
                                ^ self.perm_y[((j + dj as i64) & 255) as usize]
                                ^ self.perm_z[((k + dk as i64) & 255) as usize])
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

                let i = p.x().floor() as i64;
                let j = p.y().floor() as i64;
                let k = p.z().floor() as i64;

                let mut c = [[[Vec3::default(); 2]; 2]; 2];
                for (di, i_item) in c.iter_mut().enumerate() {
                    for (dj, ij_item) in i_item.iter_mut().enumerate() {
                        for (dk, ijk_item) in ij_item.iter_mut().enumerate() {
                            *ijk_item = self.ranvec[(self.perm_x
                                [((i + di as i64) & 255) as usize]
                                ^ self.perm_y[((j + dj as i64) & 255) as usize]
                                ^ self.perm_z[((k + dk as i64) & 255) as usize])
                                as usize]
                        }
                    }
                }

                Self::perlin_interp(c, u, v, w)
            }

            _ => {
                panic!("No matched nosie type");
            }
        }
    }

    fn perlin_generate_perm() -> Vec<i64> {
        let mut p: Vec<i64> = (0..Self::POINT_COUNT).map(|i| i as i64).collect();

        Self::permute(&mut p);

        p
    }

    fn permute(p: &mut Vec<i64>) {
        let mut rnd = rand::thread_rng();

        for i in (1..p.len()).rev() {
            let target = rnd.gen_range::<usize>(0, i);
            p.swap(i, target);
        }
    }

    fn trilinear_interp(c: [[[f64; 2]; 2]; 2], u: f64, v: f64, w: f64) -> f64 {
        let mut accum = 0.0;
        for (i, i_item) in c.iter().enumerate() {
            for (j, ij_item) in i_item.iter().enumerate() {
                for (k, ijk_item) in ij_item.iter().enumerate() {
                    accum += (i as f64 * u + (1 - i) as f64 * (1.0 - u))
                        * (j as f64 * v + (1 - j) as f64 * (1.0 - v))
                        * (k as f64 * w + (1 - k) as f64 * (1.0 - w))
                        * ijk_item;
                }
            }
        }

        accum
    }

    fn perlin_interp(c: [[[Vec3; 2]; 2]; 2], u: f64, v: f64, w: f64) -> f64 {
        let uu = u * u * (3.0 - 2.0 * u);
        let vv = v * v * (3.0 - 2.0 * v);
        let ww = w * w * (3.0 - 2.0 * w);
        let mut accum = 0.0;

        for (i, i_item) in c.iter().enumerate() {
            for (j, ij_item) in i_item.iter().enumerate() {
                for (k, ijk_item) in ij_item.iter().enumerate() {
                    let weight_v = Vec3::new(u - i as f64, v - j as f64, w - k as f64);
                    accum += (i as f64 * uu + (1.0 - i as f64) * (1.0 - uu))
                        * (j as f64 * vv + (1.0 - j as f64) * (1.0 - vv))
                        * (k as f64 * ww + (1.0 - k as f64) * (1.0 - ww))
                        * weight_v.dot(ijk_item);
                }
            }
        }

        accum
    }

    pub fn turb(&self, p: Vec3, depth: usize) -> f64 {
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
    scale: f64,
}

impl NoiseTexture {
    pub fn new(noise_type: NoiseType, scale: f64) -> Self {
        Self {
            noise: Perlin::new(noise_type),
            scale,
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, _u: f64, _v: f64, p: Vec3) -> Vec3 {
        match self.noise.noise_type {
            NoiseType::Net => Vec3::new(1.0, 1.0, 1.0) * self.noise.turb(p * self.scale, 7),

            NoiseType::Marble => {
                Vec3::new(1.0, 1.0, 1.0)
                    * 0.5
                    * (1.0 + f64::sin(self.scale * p.z() + 10.0 * self.noise.turb(p, 7)))
            }

            _ => Vec3::new(1.0, 1.0, 1.0) * 0.5 * (1.0 + self.noise.noise(p * self.scale)),
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
    fn value(&self, u: f64, v: f64, _p: Vec3) -> Vec3 {
        // If we have no texture data, then return solid cyan as a debugging aid
        if self.data.is_empty() {
            return Vec3::new(1.0, 0.0, 1.0);
        }

        // Clamp input texture coordinates to [0,1] x [1,0]
        let uu = clamp(u, 0.0, 1.0);
        let vv = 1.0 - clamp(v, 0.0, 1.0); // Filp V to image coordinates

        let mut i = (uu * f64::from(self.width)) as u32;
        let mut j = (vv * f64::from(self.height)) as u32;

        // Clamp intefer mapping, since actual coordinates should be less than 1.0
        if i >= self.width {
            i = self.width - 1;
        }
        if j >= self.height {
            j = self.height - 1;
        }

        let color_scale = 1.0 / 255.0;
        let pixel = (j * self.bytes_per_scanline + i * self.bytes_per_pixel) as usize;

        Vec3::new(
            color_scale * self.data[pixel] as f64,
            color_scale * self.data[pixel + 1] as f64,
            color_scale * self.data[pixel + 2] as f64,
        )
    }
}
