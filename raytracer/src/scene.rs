use std::sync::Arc;

use crate::color::RGB;
use crate::hittable::HittableList;

pub struct Scene {
    pub aggregate: Arc<HittableList>,
    pub lights: Arc<HittableList>,
    pub background_color: RGB,
}
