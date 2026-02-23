# cursor_playground
cursor_playground

## Generate sphere mesh (STL/3MF, solid/void/split)

This repository includes a Python script to generate watertight sphere meshes
in millimeter units.

It supports:

- Solid sphere
- Sphere with a centered spherical void
- Sphere with a centered cube void
- Split modes (`up` / `down` / `both`) after internal-void subtraction
- 2D array of generated objects (X by Y) with configurable center spacing
- Optional corner pyramid anchors at extended array corners
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
- `up`: upper half (`Z >= 0`) with split face
- `down`: lower half (`Z <= 0`) with split face
- `both`: generate both upper and lower halves

When a void is enabled (`--void-shape sphere|cube`), splitting is applied *after*
void subtraction.

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

Corner anchor option:

- `--corner-anchors`: add 4 pyramid anchors at one spacing outside the array
  bounding corners
  - Anchor footprint: `5 mm x 5 mm`
  - Anchor height: half of sphere height (`sphere_diameter / 2`)
  - Anchor bottom Z: same as sphere bottom Z

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

Split into both halves after cube void subtraction:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 20 \
  --void-shape cube \
  --void-size 8 \
  --split-half both \
  --format 3mf \
  --output sphere_20mm_void_cube_split_both.3mf
```

Array with corner anchors:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 6 \
  --array-x 2 \
  --array-y 2 \
  --array-spacing 10 \
  --corner-anchors \
  --output sphere_6mm_arr_2x2_anchors.stl
```

3MF export example:

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --format 3mf \
  --output sphere_3mm_solid.3mf
```

3MF body behavior:

- Arrays: each sphere instance is an individual 3MF body.
- `--split-half both`: each sphere contributes **two** bodies (`up` + `down`).
- `--corner-anchors`: each anchor is also its own body.

Geometry constraints:

- Sphere void diameter must be smaller than sphere diameter.
- Cube must fit inside sphere:
  `void_cube_edge * sqrt(3) < sphere_diameter`.
- If `--corner-anchors` is enabled, `--array-spacing` must be greater than `0`.

## Web UI: live preview and STL/3MF download

Open this page in a browser:

- `sphere_project/index.html`

The page allows the user to:

1. Enter sphere diameter in millimeters.
2. Choose split mode (`None`, `Upper half`, `Lower half`, `Both halves`).
3. Choose internal void shape (`None`, `Sphere`, or `Cube`).
4. Enter void size (enabled when shape is not `None` and split is `None`).
5. Set array counts (`X`, `Y`) and spacing (center-to-center, mm).
6. Optionally enable corner anchors.
7. Choose export format (`STL` or `3MF`).
8. View a live preview window.
9. Click **Generate and Download**.
10. Download the generated mesh file.

For a local web server, you can run:

```bash
cd sphere_project
python3 -m http.server 8000
```

Then open `http://localhost:8000/index.html`.
