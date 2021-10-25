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

#[derive(Clone)]
pub struct Perlin {
    ranfloat: Vec<f64>,
    perm_x: Vec<i64>,
    perm_y: Vec<i64>,
    perm_z: Vec<i64>,
}

impl Perlin {
    const POINT_COUNT: usize = 256;

    pub fn new() -> Self {
        let mut rng = rand::thread_rng();

        Self {
            ranfloat: (0..Self::POINT_COUNT)
                .map(|_| rng.gen_range::<f64>(0.0, 1.0))
                .collect(),
            perm_x: Self::perlin_generate_perm(),
            perm_y: Self::perlin_generate_perm(),
            perm_z: Self::perlin_generate_perm(),
        }
    }

    pub fn noise(&self, p: Vec3) -> f64 {
        let i: usize = ((4.0 * p.x()) as i64 & 255).try_into().unwrap();
        let j: usize = ((4.0 * p.y()) as i64 & 255).try_into().unwrap();
        let k: usize = ((4.0 * p.z()) as i64 & 255).try_into().unwrap();

        return self.ranfloat[(self.perm_x[i] ^ self.perm_y[j] ^ self.perm_z[k]) as usize];
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
}

#[derive(Clone)]
pub struct NoiseTexture {
    noise: Perlin,
}

impl NoiseTexture {
    pub fn new() -> Self {
        Self {
            noise: Perlin::new(),
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, _u: f64, _v: f64, p: Vec3) -> Vec3 {
        Vec3::new(1.0, 1.0, 1.0) * self.noise.noise(p)
    }
}
