use crate::vec3::Vec3;

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

    pub fn new_from_value(red: f64, green: f64, blue: f64) -> Self {
        Self { color_value: Vec3::new(red, green, blue) }
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
    even: SolidColor 
}

impl CheckerTexture {
    pub fn new(color1: Vec3, color2: Vec3) -> Self {
        Self { odd: SolidColor::new(color1), even: SolidColor::new(color2) }
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
