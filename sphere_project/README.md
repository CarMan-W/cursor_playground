# cursor_playground
cursor_playground

## Generate sphere STL (solid or with internal void)

This repository includes a Python script to generate a watertight STL mesh
for a sphere (surface representation) in millimeter units.

It supports:

- Solid sphere
- Sphere with a centered spherical void
- Sphere with a centered cube void
- 2D array of generated objects (X by Y) with configurable center spacing

It also includes a browser-based generator page where users can input sphere
diameter, choose internal void options, and download a generated STL file directly.

### Run

```bash
cd sphere_project
python3 generate_sphere_stl.py
```

By default this creates a solid sphere:

- `sphere_3mm_solid.stl`
- Diameter: `3.0 mm` (radius `1.5 mm`)

### Optional parameters

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --lat-segments 48 \
  --lon-segments 96 \
  --void-shape none \
  --void-size 1.0 \
  --array-x 1 \
  --array-y 1 \
  --array-spacing 100.0 \
  --output sphere_3mm_solid.stl
```

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

Geometry constraints:

- Sphere void diameter must be smaller than sphere diameter.
- Cube must fit inside sphere:
  `void_cube_edge * sqrt(3) < sphere_diameter`.

## Web UI: input dimensions and download STL

Open this page in a browser:

- `sphere_project/index.html`

The page allows the user to:

1. Enter sphere diameter in millimeters.
2. Choose internal void shape (`None`, `Sphere`, or `Cube`).
3. Enter void size (enabled when shape is not `None`).
4. Set array counts (`X`, `Y`) and spacing (center-to-center, mm).
5. Click **Generate and Download STL**.
6. Download the generated ASCII STL file.

For a local web server, you can run:

```bash
cd sphere_project
python3 -m http.server 8000
```

Then open `http://localhost:8000/index.html`.
