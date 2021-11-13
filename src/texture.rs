use std::convert::TryInto;

use crate::vec3::Vec3;
use rand::Rng;

pub trait TextureClone {
    fn clone_texture<'a>(&self) -> Box<dyn Texture>;
}

impl<T> TextureClone for T
where
    T: Texture + Clone + 'static,
{
    fn clone_texture(&self) -> Box<dyn Texture> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Texture> {
    fn clone(&self) -> Self {
        self.clone_texture()
    }
}

pub trait Texture: TextureClone {
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
pub enum NoiseType {
    // Unfiltered noise
    Square,
    // Smooth interpolation with random unit vectors
    Smooth,
    // Marble pattern with adjustable phase
    Marble,
    // Turbulent pattern that resembles a net
    Net,
    // Trillinear interpolation
    Trillinear,
}

#[derive(Clone)]
pub struct Perlin {
    ranfloat: Vec<f64>,
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

                return self.ranfloat[(self.perm_x[i] ^ self.perm_y[j] ^ self.perm_z[k]) as usize];
            }

            NoiseType::Smooth => {
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
                for di in 0..2 {
                    for dj in 0..2 {
                        for dk in 0..2 {
                            c[di][dj][dk] = self.ranfloat[(self.perm_x[(i + di as i64 & 255) as usize]
                                ^ self.perm_y[(j + dj as i64 & 255) as usize]
                                ^ self.perm_z[(k + dk as i64 & 255) as usize])
                                as usize]
                        }
                    }
                }

                Self::trilinear_interp(c, u, v, w)
            }

            _ => {
                1.0
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
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    accum += (i as f64 * u + (1 - i) as f64 * (1.0 - u))
                        * (j as f64 * v + (1 - j) as f64 * (1.0 - v))
                        * (k as f64 * w + (1 - k) as f64 * (1.0 - w))
                        * c[i][j][k];
                }
            }
        }

        accum
    }
}

#[derive(Clone)]
pub struct NoiseTexture {
    noise: Perlin,
}

impl NoiseTexture {
    pub fn new(noise_type: NoiseType) -> Self {
        Self {
            noise: Perlin::new(noise_type),
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, _u: f64, _v: f64, p: Vec3) -> Vec3 {
        Vec3::new(1.0, 1.0, 1.0) * self.noise.noise(p)
    }
}
