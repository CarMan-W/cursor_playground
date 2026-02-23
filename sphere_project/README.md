# cursor_playground
cursor_playground

## Generate sphere STL (solid or with internal void)

This repository includes a Python script to generate a watertight STL mesh
for a sphere (surface representation) in millimeter units.

It supports:

- Solid sphere
- Sphere with a centered spherical void
- Sphere with a centered cube void

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
  --output sphere_3mm_solid.stl
```

`--void-shape` options:

- `none` (default): solid sphere
- `sphere`: centered spherical void
- `cube`: centered cube void

`--void-size` meaning:

- For `sphere`: void sphere diameter (mm)
- For `cube`: void cube edge length (mm)

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
4. Click **Generate and Download STL**.
5. Download the generated ASCII STL file.

For a local web server, you can run:

```bash
cd sphere_project
python3 -m http.server 8000
```

Then open `http://localhost:8000/index.html`.
