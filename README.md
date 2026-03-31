# Yet Another Ray Tracer in Rust

![david](/output/david.png)

Ray tracer implementation in Rust, based on [Ray Tracing in One Weekend Book Series](https://raytracing.github.io/).

## Featured:
  - Spectral color [1]
  - Load and render .obj file
  - Multi-threaded Rendering
  - Static Dispatch
  - Code generation
  - QBVH with simd instructions [2]

## Usage

This binary now uses named CLI scene presets instead of a hardcoded scene selection.

```bash
cargo run --release -- --scene david
```

By default, renders are written to the scene's preset filename under `output/`.

### Available scenes

- `random-scene`
- `two-spheres`
- `two-perlin-spheres`
- `earth`
- `simple-light`
- `cornell-box`
- `cornell-box-smoke`
- `next-week-final`
- `teapot`
- `bunny`
- `three-spheres`
- `sycee`
- `david`

### Common arguments

- `--scene <name>`
- `--output <path>`
- `--width <pixels>`
- `--height <pixels>`
- `--samples <count>`
- `--max-depth <count>`
- `--workers <count>`
- `--vfov <degrees>`
- `--aperture <value>`

If only `--width` or `--height` is provided, the other dimension is derived from the selected scene's default aspect ratio.

### Examples

Render the default `david` preset:

```bash
cargo run --release -- --scene david
```

Render a smaller Cornell box image:

```bash
cargo run --release -- --scene cornell-box --width 400 --samples 32
```

Write to a custom path with camera/render overrides:

```bash
cargo run --release -- \
  --scene simple-light \
  --output renders/simple-light.png \
  --width 800 \
  --samples 128 \
  --max-depth 25 \
  --workers 8 \
  --vfov 25 \
  --aperture 0.05
```

Show the generated help text:

```bash
cargo run -- --help
```

## Bibliography
1. Smits, B. (2000). *An RGB to Spectrum Conversion for Reflectances*. University of Utah.
2. Dammertz, H., Hanika, J., & Keller, A. (2008). *Shallow Bounding Volume Hierarchies for Fast SIMD Ray Tracing of Incoherent Rays*. Computer Graphics Forum, 27(4).
