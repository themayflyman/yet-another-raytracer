use rand::Rng;
use std::cmp::PartialEq;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

#[derive(Default, Debug, Clone, Copy)]
pub struct Vec3 {
    pub elements: [f64; 3],
}

impl Vec3 {
    pub fn new(e0: f64, e1: f64, e2: f64) -> Vec3 {
        Vec3 {
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

    pub fn r(&self) -> f64 {
        self.elements[0]
    }

    pub fn g(&self) -> f64 {
        self.elements[1]
    }

    pub fn b(&self) -> f64 {
        self.elements[2]
    }

    pub fn random(low: f64, high: f64) -> Self {
        let mut rng = rand::thread_rng();
        Vec3 {
            elements: [
                rng.gen_range(low, high),
                rng.gen_range(low, high),
                rng.gen_range(low, high),
            ],
        }
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, other: Vec3) -> Self::Output {
        Vec3 {
            elements: [
                self.x() + other.x(),
                self.y() + other.y(),
                self.z() + other.z(),
            ],
        }
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, other: Vec3) -> Self::Output {
        Vec3 {
            elements: [
                self.x() - other.x(),
                self.y() - other.y(),
                self.z() - other.z(),
            ],
        }
    }
}

impl Mul<Vec3> for Vec3 {
    type Output = Self;

    fn mul(self, other: Vec3) -> Self::Output {
        Vec3 {
            elements: [
                self.x() * other.x(),
                self.y() * other.y(),
                self.z() * other.z(),
            ],
        }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Vec3 {
            elements: [self.x() * other, self.y() * other, self.z() * other],
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Self::Output {
        Vec3 {
            elements: [self * other.x(), self * other.y(), self * other.z()],
        }
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, other: f64) -> Self::Output {
        if other == 0.0 {
            return Vec3 {
                elements: [f64::MAX, f64::MAX, f64::MAX],
            };
        }
        Vec3 {
            elements: [self.x() / other, self.y() / other, self.z() / other],
        }
    }
}

impl Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Vec3 {
            elements: [-self.x(), -self.y(), -self.z()],
        }
    }
}

impl PartialEq for Vec3 {
    fn eq(&self, other: &Vec3) -> bool {
        self.x() == other.x() && self.y() == other.y() && self.z() == other.z()
    }
}

impl Index<usize> for Vec3 {
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

impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, idx: usize) -> &mut f64 {
        match idx {
            0 => &mut self.elements[0],
            1 => &mut self.elements[1],
            2 => &mut self.elements[2],
            _ => panic!("out of range"),
        }
    }
}

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

    pub fn colorize(&mut self) -> () {
        self.elements[0] = self.elements[0].floor();
        self.elements[1] = self.elements[1].floor();
        self.elements[2] = self.elements[2].floor();
    }

    pub fn near_zero(&self) -> bool {
        return f64::abs(self.elements[0]) < 1e-8
            && f64::abs(self.elements[0]) < 1e-8
            && f64::abs(self.elements[0]) < 1e-8;
    }
}
