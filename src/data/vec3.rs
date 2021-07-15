use std::cmp::PartialEq;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy)]
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
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Vec3 {
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
    type Output = Vec3;

    fn sub(self, other: Vec3) -> Vec3 {
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
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Vec3 {
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
    type Output = Vec3;

    fn mul(self, other: f64) -> Vec3 {
        Vec3 {
            elements: [self.x() * other, self.y() * other, self.z() * other],
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Vec3 {
        Vec3 {
            elements: [self * other.x(), self * other.y(), self * other.z()],
        }
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, other: f64) -> Vec3 {
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

impl PartialEq for Vec3 {
    fn eq(&self, other: &Vec3) -> bool {
        self.x() == other.x() && self.y() == other.y() && self.z() == other.z()
    }
}

impl Vec3 {
    pub fn length(&self) -> f64 {
        (self.x() * self.x() + self.y() * self.y() + self.z() * self.z()).sqrt()
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
}
