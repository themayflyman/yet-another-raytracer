// use rand::Rng;
// use std::cmp::PartialEq;
// use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

#[derive(Default, Debug, Clone, Copy)]
pub struct Vec3 {
    elements: [f64; 3],
}

#[macro_export]
macro_rules! implement_vec3_alike {
    ($type:tt) => (
        impl $type {
            pub fn new(e0: f64, e1: f64, e2: f64) -> $type {
                $type {
                    elements: [e0, e1, e2],
                }
            }

            pub fn x(&self) -> f64 {
                self.elements[0]
            }

            pub fn y(&self) -> f64 {
                self.elements[1]
            }

            pub fn z(&self) -> f64 {
                self.elements[2]
            }

            pub fn random(low: f64, high: f64) -> Self {
                use rand::Rng;
                let mut rng = rand::thread_rng();
                $type {
                    elements: [
                        rng.gen_range(low, high),
                        rng.gen_range(low, high),
                        rng.gen_range(low, high),
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
                    elements: [
                        self.x() + rhs.x(),
                        self.y() + rhs.y(),
                        self.z() + rhs.z(),
                    ],
                }
            }
        }

        impl std::ops::Sub for $type {
            type Output = Self;

            fn sub(self, rhs: $type) -> Self::Output {
                $type {
                    elements: [
                        self.x() - rhs.x(),
                        self.y() - rhs.y(),
                        self.z() - rhs.z(),
                    ],
                }
            }
        }

        impl std::ops::Mul<$type> for $type {
            type Output = Self;

            fn mul(self, rhs: $type) -> Self::Output {
                $type {
                    elements: [
                        self.x() * rhs.x(),
                        self.y() * rhs.y(),
                        self.z() * rhs.z(),
                    ],
                }
            }
        }

        impl std::ops::Mul<f64> for $type {
            type Output = Self;

            fn mul(self, other: f64) -> Self::Output {
                $type {
                    elements: [self.x() * other, self.y() * other, self.z() * other],
                }
            }
        }

        impl std::ops::Mul<$type> for f64 {
            type Output = $type;

            fn mul(self, other: $type) -> Self::Output {
                $type {
                    elements: [self * other.x(), self * other.y(), self * other.z()],
                }
            }
        }

        impl std::ops::Div<f64> for $type {
            type Output = Self;

            fn div(self, other: f64) -> Self::Output {
                if other == 0.0 {
                    return $type {
                        elements: [f64::MAX, f64::MAX, f64::MAX],
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
            type Output = f64;
            fn index(&self, idx: usize) -> &f64 {
                match idx {
                    0 => &self.elements[0],
                    1 => &self.elements[1],
                    2 => &self.elements[2],
                    _ => panic!("out of range"),
                }
            }
        }

        impl std::ops::IndexMut<usize> for $type {
            fn index_mut(&mut self, idx: usize) -> &mut f64 {
                match idx {
                    0 => &mut self.elements[0],
                    1 => &mut self.elements[1],
                    2 => &mut self.elements[2],
                    _ => panic!("out of range"),
                }
            }
        }
    )
}

implement_vec3_alike!(Vec3);

impl Vec3 {
    pub fn length(&self) -> f64 {
        (self.x() * self.x() + self.y() * self.y() + self.z() * self.z()).sqrt()
    }

    pub fn length_squared(&self) -> f64 {
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

    pub fn dot(&self, other: &Vec3) -> f64 {
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
        f64::abs(self.elements[0]) < 1e-8
            && f64::abs(self.elements[0]) < 1e-8
            && f64::abs(self.elements[0]) < 1e-8
    }
}
