mod data;
use data::vec3::Vec3;

fn main() {
    const IMAGE_WIDTH: i32 = 256;
    const IMAGE_HEIGHT: i32 = 256;

    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");

    for j in (0..IMAGE_HEIGHT).rev() {
        eprintln!("Scanlines remaining: {}", j);
        for i in 0..IMAGE_WIDTH {
            let r: f64 = i as f64 / (IMAGE_WIDTH - 1) as f64;
            let g: f64 = j as f64 / (IMAGE_HEIGHT - 1) as f64;
            let b = 0.25;

            let ir  = 255.999 * r;
            let ig = 255.999 * g;
            let ib = 255.999 * b;

            let mut color: Vec3 = Vec3::new(ir, ig, ib);
            color.colorize();

            println!("{} {} {}", color.r(), color.g(), color.b());
        }
    }

    eprintln!("Done.");
}
