# cursor_playground
cursor_playground

## Generate a 3 mm solid sphere STL

This repository includes a Python script to generate a watertight STL mesh
for a solid sphere (surface representation) in millimeter units.

### Run

```bash
python3 generate_sphere_stl.py
```

By default this creates:

- `sphere_3mm_solid.stl`
- Diameter: `3.0 mm` (radius `1.5 mm`)

### Optional parameters

```bash
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --lat-segments 48 \
  --lon-segments 96 \
  --output sphere_3mm_solid.stl
```
