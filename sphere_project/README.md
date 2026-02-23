# cursor_playground
cursor_playground

## Generate sphere mesh (STL/3MF, solid/void/split)

This repository includes a Python script to generate watertight sphere meshes
in millimeter units.

It supports:

- Solid sphere
- Sphere with a centered spherical void
- Sphere with a centered cube void
- Split-half sphere output (`up`/`down`)
- 2D array of generated objects (X by Y) with configurable center spacing
- STL or 3MF export

It also includes a browser-based generator page where users can input sphere
diameter, preview geometry, choose internal options, and download STL/3MF
directly.

### Run

```bash
cd sphere_project
python3 generate_sphere_stl.py
```

By default this creates a full solid sphere in STL:

- `sphere_3mm_solid.stl`
- Diameter: `3.0 mm` (radius `1.5 mm`)

### Optional parameters

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --lat-segments 48 \
  --lon-segments 96 \
  --split-half none \
  --void-shape none \
  --void-size 1.0 \
  --array-x 1 \
  --array-y 1 \
  --array-spacing 100.0 \
  --format auto \
  --output sphere_3mm_solid.stl
```

`--split-half` options:

- `none` (default): full sphere
- `up`: upper half sphere (`Z >= 0`) with a flat, capped split face
- `down`: lower half sphere (`Z <= 0`) with a flat, capped split face

`--void-shape` options:

- `none` (default): solid sphere
- `sphere`: centered spherical void
- `cube`: centered cube void

`--void-size` meaning:

- For `sphere`: void sphere diameter (mm)
- For `cube`: void cube edge length (mm)

Array options:

- `--array-x`: number of objects along X axis (default `1`)
- `--array-y`: number of objects along Y axis (default `1`)
- `--array-spacing`: center-to-center distance in mm for X/Y array layout
  (default `100.0`)

Output format options:

- `--format auto` (default): infer by output extension (`.3mf` => 3MF, else STL)
- `--format stl`: force ASCII STL output
- `--format 3mf`: force 3MF output

### Examples

Sphere with a centered spherical void:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --void-shape sphere \
  --void-size 1.2 \
  --output sphere_3mm_void_sphere_1p2mm.stl
```

Sphere with a centered cube void:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --void-shape cube \
  --void-size 1.0 \
  --output sphere_3mm_void_cube_1p0mm.stl
```

2D array example (10 x 10, 100 mm center spacing):

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --array-x 10 \
  --array-y 10 \
  --array-spacing 100.0 \
  --output sphere_3mm_array_10x10_100mm.stl
```

This generates `100` spheres with center-to-center distance `100 mm`
in both X and Y directions.

Split upper half example:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 20 \
  --split-half up \
  --output sphere_20mm_half_up.stl
```

3MF export example:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --format 3mf \
  --output sphere_3mm_solid.3mf
```

When exporting 3MF with an array (e.g. 10×10), each sphere is written as an
individual object/body in the file so slicers can treat them separately.

Geometry constraints:

- Sphere void diameter must be smaller than sphere diameter.
- Cube must fit inside sphere:
  `void_cube_edge * sqrt(3) < sphere_diameter`.
- When `--split-half` is `up` or `down`, internal void must be `none`.

## Web UI: live preview and STL/3MF download

Open this page in a browser:

- `sphere_project/index.html`

The page allows the user to:

1. Enter sphere diameter in millimeters.
2. Choose split mode (`None`, `Upper half`, `Lower half`).
3. Choose internal void shape (`None`, `Sphere`, or `Cube`).
4. Enter void size (enabled when shape is not `None` and split is `None`).
5. Set array counts (`X`, `Y`) and spacing (center-to-center, mm).
6. Choose export format (`STL` or `3MF`).
7. View a live preview window.
8. Click **Generate and Download**.
9. Download the generated mesh file.

For a local web server, you can run:

```bash
cd sphere_project
python3 -m http.server 8000
```

Then open `http://localhost:8000/index.html`.
