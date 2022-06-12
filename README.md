# Yet Another Ray Tracer in Rust

Ray tracer implementation in Rust, based on [Ray Tracing in One Weekend Book Series](https://raytracing.github.io/).

## Featured:
  - Multi-threaded Rendering
  - Static Dispatch
  - Code generation

which siginificantly speed up the rendering, the following picture took 30 workers 2635 seconds to render on my 5950X workstation:
![the next week final scene](/output/the_next_week_final_static_scene.png)

## Usages
Modify scene no. in `main.rs` or implement your own scene in `scenes.rs`. Then run
```
cargo run --release
```

## Showcases

![book1 final scene](/output/random_scene.png)

![cornel box that uses mixture pdf to reduce noise](/output/cornell_box.png)

![linear gradient texture](/output/two_perlin_spheres.png)
