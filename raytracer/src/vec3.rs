use core::f32;

// use rand::Rng;
// use std::cmp::PartialEq;
// use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use rand::Rng;

#[derive(Default, Debug, Clone, Copy)]
pub struct Vec3 {
    elements: [f32; 3],
}

#[macro_export]
macro_rules! implement_vec3_alike {
    ($type:tt) => {
        impl $type {
            pub fn new(e0: f32, e1: f32, e2: f32) -> $type {
                $type {
                    elements: [e0, e1, e2],
                }
            }

            pub fn x(&self) -> f32 {
                self.elements[0]
            }

            pub fn y(&self) -> f32 {
                self.elements[1]
            }

            pub fn z(&self) -> f32 {
                self.elements[2]
            }

            pub fn random(low: f32, high: f32) -> Self {
                use rand::Rng;
                let mut rng = rand::thread_rng();
                $type {
                    elements: [
                        rng.gen_range(low..high),
                        rng.gen_range(low..high),
                        rng.gen_range(low..high),
                    ],
                }
            }
        }

        impl std::ops::Add for $type {
            type Output = Self;

            fn add(self, other: $type) -> Self::Output {
                $type {
                    elements: [
                        self.x() + other.x(),
                        self.y() + other.y(),
                        self.z() + other.z(),
                    ],
                }
            }
        }

        impl std::ops::AddAssign for $type {
            fn add_assign(&mut self, rhs: Self) {
                *self = Self {
                    elements: [self.x() + rhs.x(), self.y() + rhs.y(), self.z() + rhs.z()],
                }
            }
        }

        impl std::ops::Sub for $type {
            type Output = Self;

            fn sub(self, rhs: $type) -> Self::Output {
                $type {
                    elements: [self.x() - rhs.x(), self.y() - rhs.y(), self.z() - rhs.z()],
                }
            }
        }

        impl std::ops::Mul<$type> for $type {
            type Output = Self;

            fn mul(self, rhs: $type) -> Self::Output {
                $type {
                    elements: [self.x() * rhs.x(), self.y() * rhs.y(), self.z() * rhs.z()],
                }
            }
        }

        impl std::ops::Mul<f32> for $type {
            type Output = Self;

            fn mul(self, other: f32) -> Self::Output {
                $type {
                    elements: [self.x() * other, self.y() * other, self.z() * other],
                }
            }
        }

        impl std::ops::Mul<$type> for f32 {
            type Output = $type;

            fn mul(self, other: $type) -> Self::Output {
                $type {
                    elements: [self * other.x(), self * other.y(), self * other.z()],
                }
            }
        }

        impl std::ops::Div<f32> for $type {
            type Output = Self;

            fn div(self, other: f32) -> Self::Output {
                if other == 0.0 {
                    return $type {
                        elements: [f32::MAX, f32::MAX, f32::MAX],
                    };
                }
                $type {
                    elements: [self.x() / other, self.y() / other, self.z() / other],
                }
            }
        }

        impl std::ops::Neg for $type {
            type Output = Self;

            fn neg(self) -> Self::Output {
                $type {
                    elements: [-self.x(), -self.y(), -self.z()],
                }
            }
        }

        impl std::cmp::PartialEq for $type {
            fn eq(&self, other: &$type) -> bool {
                self.x() == other.x() && self.y() == other.y() && self.z() == other.z()
            }
        }

        impl std::ops::Index<usize> for $type {
            type Output = f32;
            fn index(&self, idx: usize) -> &f32 {
                match idx {
                    0 => &self.elements[0],
                    1 => &self.elements[1],
                    2 => &self.elements[2],
                    _ => panic!("out of range"),
                }
            }
        }

        impl std::ops::IndexMut<usize> for $type {
            fn index_mut(&mut self, idx: usize) -> &mut f32 {
                match idx {
                    0 => &mut self.elements[0],
                    1 => &mut self.elements[1],
                    2 => &mut self.elements[2],
                    _ => panic!("out of range"),
                }
            }
        }

        impl std::convert::From<std::simd::f32x4> for $type {
            fn from(v: std::simd::f32x4) -> Self {
                $type::new(v[0], v[1], v[2])
            }
        }

        impl std::convert::From<$type> for std::simd::f32x4 {
            fn from(v: $type) -> Self {
                std::simd::f32x4::from_array([v.x(), v.y(), v.z(), v.z()])
            }
        }

        impl std::convert::From<$type> for std::simd::f32x8 {
            fn from(v: $type) -> Self {
                std::simd::f32x8::from_array([
                    v.x(),
                    v.y(),
                    v.z(),
                    v.z(),
                    v.x(),
                    v.y(),
                    v.z(),
                    v.z(),
                ])
            }
        }
    };
}

implement_vec3_alike!(Vec3);

impl Vec3 {
    pub fn length(&self) -> f32 {
        (self.x() * self.x() + self.y() * self.y() + self.z() * self.z()).sqrt()
    }

    pub fn length_squared(&self) -> f32 {
        self.x() * self.x() + self.y() * self.y() + self.z() * self.z()
    }

    pub fn unit_vector(&self) -> Vec3 {
        Vec3 {
            elements: [
                self.x() / self.length(),
                self.y() / self.length(),
                self.z() / self.length(),
            ],
        }
    }

    pub fn random_unit_vector<R: Rng>(rng: &mut R) -> Self {
        let a = rng.gen_range(0.0..2.0 * f32::consts::PI);
        let z = rng.gen_range(-1.0..1.0);
        let r = f32::sqrt(1.0 - z * z);

        Vec3::new(r * a.cos(), r * a.sin(), z)
    }

    pub fn dot(&self, other: &Vec3) -> f32 {
        (self.x() * other.x()) + (self.y() * other.y()) + (self.z() * other.z())
    }

    pub fn cross(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            elements: [
                self.y() * other.z() - self.z() * other.y(),
                self.z() * other.x() - self.x() * other.z(),
                self.x() * other.y() - self.y() * other.x(),
            ],
        }
    }

    pub fn near_zero(&self) -> bool {
        f32::abs(self.elements[0]) < 1e-8
            && f32::abs(self.elements[0]) < 1e-8
            && f32::abs(self.elements[0]) < 1e-8
    }

    pub fn abs(&self) -> Vec3 {
        Vec3 {
            elements: [self.x().abs(), self.y().abs(), self.z().abs()],
        }
    }
}
