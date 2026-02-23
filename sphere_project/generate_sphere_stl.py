#!/usr/bin/env python3
"""Generate watertight sphere meshes in millimeters.

Supports:
- Solid sphere
- Sphere with centered internal void (sphere/cube)
- Split sphere modes (up/down/both), including split-after-void
- 2D X/Y array layout with configurable center spacing
- Optional corner pyramid anchors
- STL or 3MF export
"""

from __future__ import annotations

import argparse
import math
import os
import zipfile
from dataclasses import dataclass
from html import escape
from typing import List, Sequence, Tuple

Vec3 = Tuple[float, float, float]


@dataclass(frozen=True)
class Triangle:
    v1: Vec3
    v2: Vec3
    v3: Vec3


def subtract(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def normalize(v: Vec3) -> Vec3:
    length = math.sqrt(dot(v, v))
    if length == 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / length, v[1] / length, v[2] / length)


def average(*points: Vec3) -> Vec3:
    count = float(len(points))
    return (
        sum(p[0] for p in points) / count,
        sum(p[1] for p in points) / count,
        sum(p[2] for p in points) / count,
    )


def orient_outward(v1: Vec3, v2: Vec3, v3: Vec3) -> Triangle:
    """Ensure triangle winding produces outward normal for origin-centered sphere."""
    n = cross(subtract(v2, v1), subtract(v3, v1))
    c = average(v1, v2, v3)
    if dot(n, c) < 0.0:
        return Triangle(v1=v1, v2=v3, v3=v2)
    return Triangle(v1=v1, v2=v2, v3=v3)


def orient_inward(v1: Vec3, v2: Vec3, v3: Vec3) -> Triangle:
    """Ensure triangle winding produces inward normal for origin-centered meshes."""
    n = cross(subtract(v2, v1), subtract(v3, v1))
    c = average(v1, v2, v3)
    if dot(n, c) > 0.0:
        return Triangle(v1=v1, v2=v3, v3=v2)
    return Triangle(v1=v1, v2=v2, v3=v3)


def build_sphere_mesh(radius_mm: float, lat_segments: int, lon_segments: int) -> List[Triangle]:
    if radius_mm <= 0.0:
        raise ValueError("Radius must be > 0.")
    if lat_segments < 2:
        raise ValueError("lat_segments must be >= 2.")
    if lon_segments < 3:
        raise ValueError("lon_segments must be >= 3.")

    top: Vec3 = (0.0, 0.0, radius_mm)
    bottom: Vec3 = (0.0, 0.0, -radius_mm)

    rings: List[List[Vec3]] = []
    for i in range(1, lat_segments):
        theta = math.pi * i / lat_segments
        z = radius_mm * math.cos(theta)
        ring_radius = radius_mm * math.sin(theta)
        ring: List[Vec3] = []
        for j in range(lon_segments):
            phi = 2.0 * math.pi * j / lon_segments
            x = ring_radius * math.cos(phi)
            y = ring_radius * math.sin(phi)
            ring.append((x, y, z))
        rings.append(ring)

    triangles: List[Triangle] = []
    first_ring = rings[0]
    last_ring = rings[-1]

    # Top cap.
    for j in range(lon_segments):
        v2 = first_ring[(j + 1) % lon_segments]
        v3 = first_ring[j]
        triangles.append(orient_outward(top, v2, v3))

    # Middle strips.
    for ring_idx in range(len(rings) - 1):
        ring_a = rings[ring_idx]
        ring_b = rings[ring_idx + 1]
        for j in range(lon_segments):
            v00 = ring_a[j]
            v01 = ring_a[(j + 1) % lon_segments]
            v10 = ring_b[j]
            v11 = ring_b[(j + 1) % lon_segments]
            triangles.append(orient_outward(v00, v01, v11))
            triangles.append(orient_outward(v00, v11, v10))

    # Bottom cap.
    for j in range(lon_segments):
        v2 = last_ring[j]
        v3 = last_ring[(j + 1) % lon_segments]
        triangles.append(orient_outward(bottom, v2, v3))

    return triangles


def build_hemisphere_mesh(
    radius_mm: float,
    lat_segments: int,
    lon_segments: int,
    split_half: str,
    include_equator_cap: bool = True,
) -> List[Triangle]:
    """Build a watertight upper or lower half sphere with an equatorial cap."""
    if radius_mm <= 0.0:
        raise ValueError("Radius must be > 0.")
    if lat_segments < 2:
        raise ValueError("lat_segments must be >= 2.")
    if lon_segments < 3:
        raise ValueError("lon_segments must be >= 3.")
    if split_half not in ("up", "down"):
        raise ValueError("split_half must be either 'up' or 'down'.")

    hemi_lat_segments = max(1, lat_segments // 2)
    pole: Vec3 = (0.0, 0.0, radius_mm if split_half == "up" else -radius_mm)

    rings: List[List[Vec3]] = []
    for i in range(1, hemi_lat_segments + 1):
        theta = (math.pi / 2.0) * i / hemi_lat_segments
        z = radius_mm * math.cos(theta)
        if split_half == "down":
            z = -z
        ring_radius = radius_mm * math.sin(theta)
        ring: List[Vec3] = []
        for j in range(lon_segments):
            phi = 2.0 * math.pi * j / lon_segments
            x = ring_radius * math.cos(phi)
            y = ring_radius * math.sin(phi)
            ring.append((x, y, z))
        rings.append(ring)

    triangles: List[Triangle] = []
    first_ring = rings[0]

    # Pole cap.
    for j in range(lon_segments):
        if split_half == "up":
            v2 = first_ring[(j + 1) % lon_segments]
            v3 = first_ring[j]
        else:
            v2 = first_ring[j]
            v3 = first_ring[(j + 1) % lon_segments]
        triangles.append(orient_outward(pole, v2, v3))

    # Middle strips.
    for ring_idx in range(len(rings) - 1):
        ring_a = rings[ring_idx]
        ring_b = rings[ring_idx + 1]
        for j in range(lon_segments):
            v00 = ring_a[j]
            v01 = ring_a[(j + 1) % lon_segments]
            v10 = ring_b[j]
            v11 = ring_b[(j + 1) % lon_segments]
            triangles.append(orient_outward(v00, v01, v11))
            triangles.append(orient_outward(v00, v11, v10))

    if include_equator_cap:
        # Equator cap keeps the half sphere watertight.
        center: Vec3 = (0.0, 0.0, 0.0)
        equator_ring = rings[-1]
        for j in range(lon_segments):
            if split_half == "up":
                triangles.append(
                    Triangle(
                        v1=center,
                        v2=equator_ring[(j + 1) % lon_segments],
                        v3=equator_ring[j],
                    )
                )
            else:
                triangles.append(
                    Triangle(
                        v1=center,
                        v2=equator_ring[j],
                        v3=equator_ring[(j + 1) % lon_segments],
                    )
                )

    return triangles


def build_cube_cavity_mesh(edge_mm: float) -> List[Triangle]:
    if edge_mm <= 0.0:
        raise ValueError("Cube void size must be > 0.")

    half = edge_mm / 2.0
    faces: Sequence[Tuple[Vec3, Vec3, Vec3, Vec3]] = (
        # +X
        ((half, -half, -half), (half, half, -half), (half, half, half), (half, -half, half)),
        # -X
        ((-half, -half, -half), (-half, -half, half), (-half, half, half), (-half, half, -half)),
        # +Y
        ((-half, half, -half), (-half, half, half), (half, half, half), (half, half, -half)),
        # -Y
        ((-half, -half, -half), (half, -half, -half), (half, -half, half), (-half, -half, half)),
        # +Z
        ((-half, -half, half), (half, -half, half), (half, half, half), (-half, half, half)),
        # -Z
        ((-half, -half, -half), (-half, half, -half), (half, half, -half), (half, -half, -half)),
    )

    triangles: List[Triangle] = []
    for a, b, c, d in faces:
        triangles.append(orient_inward(a, b, c))
        triangles.append(orient_inward(a, c, d))
    return triangles


def validate_internal_void(outer_radius_mm: float, void_shape: str, void_size_mm: float) -> None:
    if void_shape == "none":
        return

    if void_size_mm <= 0.0:
        raise ValueError("void_size must be > 0 when using an internal void.")

    if void_shape == "sphere":
        void_radius = void_size_mm / 2.0
        if void_radius >= outer_radius_mm:
            raise ValueError(
                "Sphere void diameter must be smaller than the outer sphere diameter."
            )
        return

    if void_shape == "cube":
        half_diagonal = (void_size_mm * math.sqrt(3.0)) / 2.0
        if half_diagonal >= outer_radius_mm:
            raise ValueError(
                "Cube void is too large; cube must fit fully inside the outer sphere."
            )
        return

    raise ValueError(f"Unsupported void_shape: {void_shape}")


def validate_split_config(split_half: str) -> None:
    if split_half not in ("none", "up", "down", "both"):
        raise ValueError("split_half must be one of: none, up, down, both.")


def orient_triangle_to_z(v1: Vec3, v2: Vec3, v3: Vec3, positive_z: bool) -> Triangle:
    n = cross(subtract(v2, v1), subtract(v3, v1))
    nz = n[2]
    if positive_z:
        if nz < 0.0:
            return Triangle(v1=v1, v2=v3, v3=v2)
    elif nz > 0.0:
        return Triangle(v1=v1, v2=v3, v3=v2)
    return Triangle(v1=v1, v2=v2, v3=v3)


def build_equator_ring(radius_mm: float, lon_segments: int) -> List[Vec3]:
    ring: List[Vec3] = []
    for j in range(lon_segments):
        phi = 2.0 * math.pi * j / lon_segments
        ring.append((radius_mm * math.cos(phi), radius_mm * math.sin(phi), 0.0))
    return ring


def build_square_boundary_ring(half_edge_mm: float, lon_segments: int) -> List[Vec3]:
    ring: List[Vec3] = []
    for j in range(lon_segments):
        phi = 2.0 * math.pi * j / lon_segments
        c = math.cos(phi)
        s = math.sin(phi)
        scale = half_edge_mm / max(abs(c), abs(s))
        ring.append((scale * c, scale * s, 0.0))
    return ring


def build_split_cap_between_rings(
    outer_ring: Sequence[Vec3],
    inner_ring: Sequence[Vec3],
    split_half: str,
) -> List[Triangle]:
    if len(outer_ring) != len(inner_ring):
        raise ValueError("Outer and inner cap rings must have the same segment count.")
    if split_half not in ("up", "down"):
        raise ValueError("split_half must be 'up' or 'down'.")

    triangles: List[Triangle] = []
    positive_z = split_half == "down"
    ring_size = len(outer_ring)
    for j in range(ring_size):
        o0 = outer_ring[j]
        o1 = outer_ring[(j + 1) % ring_size]
        i0 = inner_ring[j]
        i1 = inner_ring[(j + 1) % ring_size]
        triangles.append(orient_triangle_to_z(o0, o1, i1, positive_z=positive_z))
        triangles.append(orient_triangle_to_z(o0, i1, i0, positive_z=positive_z))
    return triangles


def build_cube_cavity_half_open_mesh(
    edge_mm: float,
    lon_segments: int,
    split_half: str,
) -> Tuple[List[Triangle], List[Vec3]]:
    if edge_mm <= 0.0:
        raise ValueError("Cube void size must be > 0.")
    if split_half not in ("up", "down"):
        raise ValueError("split_half must be 'up' or 'down'.")

    half = edge_mm / 2.0
    split_ring = build_square_boundary_ring(half, lon_segments)
    z_far = half if split_half == "up" else -half
    far_ring: List[Vec3] = [(x, y, z_far) for x, y, _ in split_ring]

    triangles: List[Triangle] = []

    # Side walls from split plane to far cube face.
    for j in range(lon_segments):
        s0 = split_ring[j]
        s1 = split_ring[(j + 1) % lon_segments]
        f0 = far_ring[j]
        f1 = far_ring[(j + 1) % lon_segments]
        triangles.append(orient_inward(s0, s1, f1))
        triangles.append(orient_inward(s0, f1, f0))

    # Far face of the void cavity.
    center_far: Vec3 = (0.0, 0.0, z_far)
    for j in range(lon_segments):
        f0 = far_ring[j]
        f1 = far_ring[(j + 1) % lon_segments]
        triangles.append(orient_inward(center_far, f0, f1))

    return triangles, split_ring


def build_split_void_half_mesh(
    radius_mm: float,
    lat_segments: int,
    lon_segments: int,
    split_half: str,
    void_shape: str,
    void_size_mm: float,
) -> List[Triangle]:
    if split_half not in ("up", "down"):
        raise ValueError("split_half must be 'up' or 'down'.")
    if void_shape == "none":
        return build_hemisphere_mesh(
            radius_mm=radius_mm,
            lat_segments=lat_segments,
            lon_segments=lon_segments,
            split_half=split_half,
        )

    validate_internal_void(radius_mm, void_shape, void_size_mm)

    outer_surface = build_hemisphere_mesh(
        radius_mm=radius_mm,
        lat_segments=lat_segments,
        lon_segments=lon_segments,
        split_half=split_half,
        include_equator_cap=False,
    )
    outer_ring = build_equator_ring(radius_mm, lon_segments)

    if void_shape == "sphere":
        inner_radius = void_size_mm / 2.0
        inner_surface = build_hemisphere_mesh(
            radius_mm=inner_radius,
            lat_segments=lat_segments,
            lon_segments=lon_segments,
            split_half=split_half,
            include_equator_cap=False,
        )
        inner_surface = [
            Triangle(v1=tri.v1, v2=tri.v3, v3=tri.v2) for tri in inner_surface
        ]
        inner_ring = build_equator_ring(inner_radius, lon_segments)
    elif void_shape == "cube":
        inner_surface, inner_ring = build_cube_cavity_half_open_mesh(
            edge_mm=void_size_mm,
            lon_segments=lon_segments,
            split_half=split_half,
        )
    else:
        raise ValueError(f"Unsupported void_shape: {void_shape}")

    split_cap = build_split_cap_between_rings(outer_ring, inner_ring, split_half)
    return outer_surface + inner_surface + split_cap


def build_single_sphere_mesh(
    radius_mm: float,
    lat_segments: int,
    lon_segments: int,
    void_shape: str,
    void_size_mm: float,
    split_half: str,
) -> List[Triangle]:
    if split_half in ("up", "down"):
        return build_split_void_half_mesh(
            radius_mm=radius_mm,
            lat_segments=lat_segments,
            lon_segments=lon_segments,
            split_half=split_half,
            void_shape=void_shape,
            void_size_mm=void_size_mm,
        )

    triangles = build_sphere_mesh(
        radius_mm=radius_mm,
        lat_segments=lat_segments,
        lon_segments=lon_segments,
    )
    validate_internal_void(radius_mm, void_shape, void_size_mm)

    if void_shape == "sphere":
        cavity = build_sphere_mesh(
            radius_mm=void_size_mm / 2.0,
            lat_segments=lat_segments,
            lon_segments=lon_segments,
        )
        for tri in cavity:
            triangles.append(Triangle(v1=tri.v1, v2=tri.v3, v3=tri.v2))
    elif void_shape == "cube":
        triangles.extend(build_cube_cavity_mesh(void_size_mm))
    elif void_shape != "none":
        raise ValueError(f"Unsupported void_shape: {void_shape}")

    return triangles


def build_sphere_body_meshes(
    radius_mm: float,
    lat_segments: int,
    lon_segments: int,
    void_shape: str,
    void_size_mm: float,
    split_half: str,
) -> List[List[Triangle]]:
    validate_split_config(split_half)
    if split_half == "both":
        return [
            build_single_sphere_mesh(
                radius_mm=radius_mm,
                lat_segments=lat_segments,
                lon_segments=lon_segments,
                void_shape=void_shape,
                void_size_mm=void_size_mm,
                split_half="up",
            ),
            build_single_sphere_mesh(
                radius_mm=radius_mm,
                lat_segments=lat_segments,
                lon_segments=lon_segments,
                void_shape=void_shape,
                void_size_mm=void_size_mm,
                split_half="down",
            ),
        ]

    return [
        build_single_sphere_mesh(
            radius_mm=radius_mm,
            lat_segments=lat_segments,
            lon_segments=lon_segments,
            void_shape=void_shape,
            void_size_mm=void_size_mm,
            split_half=split_half,
        )
    ]


def validate_array_config(array_x: int, array_y: int, array_spacing_mm: float) -> None:
    if array_x < 1:
        raise ValueError("array_x must be >= 1.")
    if array_y < 1:
        raise ValueError("array_y must be >= 1.")
    if (array_x > 1 or array_y > 1) and array_spacing_mm <= 0.0:
        raise ValueError("array_spacing must be > 0 when array_x or array_y is greater than 1.")


def offset_point(v: Vec3, dx: float, dy: float, dz: float) -> Vec3:
    return (v[0] + dx, v[1] + dy, v[2] + dz)


def offset_triangle(tri: Triangle, dx: float, dy: float, dz: float) -> Triangle:
    return Triangle(
        v1=offset_point(tri.v1, dx, dy, dz),
        v2=offset_point(tri.v2, dx, dy, dz),
        v3=offset_point(tri.v3, dx, dy, dz),
    )


def tile_mesh_2d(
    base_triangles: Sequence[Triangle],
    array_x: int,
    array_y: int,
    array_spacing_mm: float,
) -> List[Triangle]:
    validate_array_config(array_x, array_y, array_spacing_mm)

    if array_x == 1 and array_y == 1:
        return list(base_triangles)

    triangles: List[Triangle] = []
    x_center = (array_x - 1) / 2.0
    y_center = (array_y - 1) / 2.0

    for yi in range(array_y):
        dy = (yi - y_center) * array_spacing_mm
        for xi in range(array_x):
            dx = (xi - x_center) * array_spacing_mm
            for tri in base_triangles:
                triangles.append(offset_triangle(tri, dx, dy, 0.0))

    return triangles


def tile_body_meshes_2d(
    base_meshes: Sequence[Sequence[Triangle]],
    array_x: int,
    array_y: int,
    array_spacing_mm: float,
) -> List[List[Triangle]]:
    validate_array_config(array_x, array_y, array_spacing_mm)

    x_center = (array_x - 1) / 2.0
    y_center = (array_y - 1) / 2.0
    tiled: List[List[Triangle]] = []
    for yi in range(array_y):
        dy = (yi - y_center) * array_spacing_mm
        for xi in range(array_x):
            dx = (xi - x_center) * array_spacing_mm
            for mesh in base_meshes:
                tiled.append([offset_triangle(tri, dx, dy, 0.0) for tri in mesh])
    return tiled


def flatten_meshes(mesh_list: Sequence[Sequence[Triangle]]) -> List[Triangle]:
    triangles: List[Triangle] = []
    for mesh in mesh_list:
        triangles.extend(mesh)
    return triangles


def orient_triangle_away_from_point(
    v1: Vec3,
    v2: Vec3,
    v3: Vec3,
    inside_point: Vec3,
) -> Triangle:
    n = cross(subtract(v2, v1), subtract(v3, v1))
    c = average(v1, v2, v3)
    to_face = subtract(c, inside_point)
    if dot(n, to_face) < 0.0:
        return Triangle(v1=v1, v2=v3, v3=v2)
    return Triangle(v1=v1, v2=v2, v3=v3)


def build_anchor_pyramid_mesh(
    center_x_mm: float,
    center_y_mm: float,
    base_z_mm: float,
    sphere_diameter_mm: float,
) -> List[Triangle]:
    half_footprint = 2.5  # 5 mm x 5 mm base footprint.
    height = sphere_diameter_mm / 2.0
    if height <= 0.0:
        raise ValueError("sphere diameter must be > 0 for anchor generation.")

    p00: Vec3 = (center_x_mm - half_footprint, center_y_mm - half_footprint, base_z_mm)
    p10: Vec3 = (center_x_mm + half_footprint, center_y_mm - half_footprint, base_z_mm)
    p11: Vec3 = (center_x_mm + half_footprint, center_y_mm + half_footprint, base_z_mm)
    p01: Vec3 = (center_x_mm - half_footprint, center_y_mm + half_footprint, base_z_mm)
    apex: Vec3 = (center_x_mm, center_y_mm, base_z_mm + height)

    inside_point: Vec3 = (center_x_mm, center_y_mm, base_z_mm + (height / 4.0))
    triangles = [
        # Side faces.
        orient_triangle_away_from_point(p00, p10, apex, inside_point),
        orient_triangle_away_from_point(p10, p11, apex, inside_point),
        orient_triangle_away_from_point(p11, p01, apex, inside_point),
        orient_triangle_away_from_point(p01, p00, apex, inside_point),
        # Base face.
        orient_triangle_away_from_point(p00, p11, p10, inside_point),
        orient_triangle_away_from_point(p00, p01, p11, inside_point),
    ]
    return triangles


def corner_anchor_positions(
    array_x: int,
    array_y: int,
    array_spacing_mm: float,
) -> List[Tuple[float, float]]:
    validate_array_config(array_x, array_y, array_spacing_mm)
    if array_spacing_mm <= 0.0:
        raise ValueError("array_spacing must be > 0 when using corner anchors.")

    x_min = -((array_x - 1) / 2.0) * array_spacing_mm
    x_max = ((array_x - 1) / 2.0) * array_spacing_mm
    y_min = -((array_y - 1) / 2.0) * array_spacing_mm
    y_max = ((array_y - 1) / 2.0) * array_spacing_mm

    return [
        (x_min - array_spacing_mm, y_min - array_spacing_mm),
        (x_max + array_spacing_mm, y_max + array_spacing_mm),
        (x_min - array_spacing_mm, y_max + array_spacing_mm),
        (x_max + array_spacing_mm, y_min - array_spacing_mm),
    ]


def build_corner_anchor_meshes(
    enabled: bool,
    sphere_diameter_mm: float,
    sphere_radius_mm: float,
    array_x: int,
    array_y: int,
    array_spacing_mm: float,
) -> List[List[Triangle]]:
    if not enabled:
        return []

    anchors: List[List[Triangle]] = []
    for cx, cy in corner_anchor_positions(array_x, array_y, array_spacing_mm):
        anchors.append(
            build_anchor_pyramid_mesh(
                center_x_mm=cx,
                center_y_mm=cy,
                base_z_mm=-sphere_radius_mm,
                sphere_diameter_mm=sphere_diameter_mm,
            )
        )
    return anchors


def triangle_normal(tri: Triangle) -> Vec3:
    raw_normal = cross(subtract(tri.v2, tri.v1), subtract(tri.v3, tri.v1))
    return normalize(raw_normal)


def write_ascii_stl(path: str, solid_name: str, triangles: Sequence[Triangle]) -> None:
    with open(path, "w", encoding="ascii") as f:
        f.write(f"solid {solid_name}\n")
        for tri in triangles:
            n = triangle_normal(tri)
            f.write(f"  facet normal {n[0]:.8f} {n[1]:.8f} {n[2]:.8f}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {tri.v1[0]:.8f} {tri.v1[1]:.8f} {tri.v1[2]:.8f}\n")
            f.write(f"      vertex {tri.v2[0]:.8f} {tri.v2[1]:.8f} {tri.v2[2]:.8f}\n")
            f.write(f"      vertex {tri.v3[0]:.8f} {tri.v3[1]:.8f} {tri.v3[2]:.8f}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")


def build_indexed_mesh(
    triangles: Sequence[Triangle],
) -> Tuple[List[Vec3], List[Tuple[int, int, int]]]:
    vertices: List[Vec3] = []
    vertex_to_index: dict[Vec3, int] = {}
    indexed_triangles: List[Tuple[int, int, int]] = []

    for tri in triangles:
        tri_indices: List[int] = []
        for vertex in (tri.v1, tri.v2, tri.v3):
            idx = vertex_to_index.get(vertex)
            if idx is None:
                idx = len(vertices)
                vertices.append(vertex)
                vertex_to_index[vertex] = idx
            tri_indices.append(idx)
        indexed_triangles.append((tri_indices[0], tri_indices[1], tri_indices[2]))

    return vertices, indexed_triangles


def write_3mf(
    path: str,
    model_name: str,
    mesh_list: Sequence[Sequence[Triangle]],
) -> None:
    """Write 3MF with one object/body per mesh in mesh_list."""
    model_lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<model unit="millimeter" xml:lang="en-US" '
        'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">',
        "  <resources>",
    ]
    for obj_id, triangles in enumerate(mesh_list, start=1):
        vertices, indexed_triangles = build_indexed_mesh(triangles)
        body_name = model_name if len(mesh_list) == 1 else f"{model_name}_{obj_id}"
        safe_name = escape(body_name, quote=True)
        model_lines.append(f'    <object id="{obj_id}" type="model" name="{safe_name}">')
        model_lines.append("      <mesh>")
        model_lines.append("        <vertices>")
        for x, y, z in vertices:
            model_lines.append(f'          <vertex x="{x:.8f}" y="{y:.8f}" z="{z:.8f}"/>')
        model_lines.append("        </vertices>")
        model_lines.append("        <triangles>")
        for i1, i2, i3 in indexed_triangles:
            model_lines.append(f'          <triangle v1="{i1}" v2="{i2}" v3="{i3}"/>')
        model_lines.append("        </triangles>")
        model_lines.append("      </mesh>")
        model_lines.append("    </object>")
    model_lines.append("  </resources>")
    model_lines.append("  <build>")
    for obj_id in range(1, len(mesh_list) + 1):
        model_lines.append(f'    <item objectid="{obj_id}"/>')
    model_lines.extend(["  </build>", "</model>", ""])
    model_xml = "\n".join(model_lines)

    content_types = """<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>
</Types>
"""

    root_rels = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>
</Relationships>
"""

    with zipfile.ZipFile(path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("3D/3dmodel.model", model_xml)


def resolve_output_format(requested_format: str, output_path: str) -> str:
    if requested_format in ("stl", "3mf"):
        return requested_format
    ext = os.path.splitext(output_path)[1].lower()
    if ext == ".3mf":
        return "3mf"
    return "stl"


def write_mesh_output(
    path: str,
    output_format: str,
    solid_name: str,
    triangles: Sequence[Triangle],
    mesh_list: Sequence[Sequence[Triangle]] | None = None,
) -> None:
    if output_format == "stl":
        write_ascii_stl(path, solid_name, triangles)
        return
    if output_format == "3mf":
        bodies = mesh_list if mesh_list is not None else [triangles]
        write_3mf(path, solid_name, bodies)
        return
    raise ValueError(f"Unsupported output format: {output_format}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate watertight sphere meshes in STL or 3MF format (millimeters)."
    )
    parser.add_argument(
        "--diameter",
        type=float,
        default=3.0,
        help="Sphere diameter in mm (default: 3.0).",
    )
    parser.add_argument(
        "--lat-segments",
        type=int,
        default=48,
        help="Number of latitude segments from top to bottom (default: 48).",
    )
    parser.add_argument(
        "--lon-segments",
        type=int,
        default=96,
        help="Number of longitude segments around axis (default: 96).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="sphere_3mm_solid.stl",
        help="Output mesh file path (.stl or .3mf).",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=("auto", "stl", "3mf"),
        default="auto",
        help=(
            "Output format. Use auto to infer from output extension "
            "(default: auto)."
        ),
    )
    parser.add_argument(
        "--void-shape",
        type=str,
        choices=("none", "sphere", "cube"),
        default="none",
        help="Internal void shape: none, sphere, or cube (default: none).",
    )
    parser.add_argument(
        "--void-size",
        type=float,
        default=1.0,
        help=(
            "Internal void size in mm. Sphere uses diameter; cube uses edge length "
            "(default: 1.0). Ignored when --void-shape=none."
        ),
    )
    parser.add_argument(
        "--array-x",
        type=int,
        default=1,
        help="Number of objects along X axis (default: 1).",
    )
    parser.add_argument(
        "--array-y",
        type=int,
        default=1,
        help="Number of objects along Y axis (default: 1).",
    )
    parser.add_argument(
        "--array-spacing",
        type=float,
        default=100.0,
        help=(
            "Center-to-center spacing in mm used for the 2D array "
            "(default: 100.0)."
        ),
    )
    parser.add_argument(
        "--corner-anchors",
        action="store_true",
        help=(
            "Add four corner pyramid anchors at one spacing outside the array bounds. "
            "Anchors are 5x5 mm footprint, half sphere height, and share the "
            "same bottom Z as the spheres."
        ),
    )
    parser.add_argument(
        "--split-half",
        type=str,
        choices=("none", "up", "down", "both"),
        default="none",
        help=(
            "Split mode: none (full), up, down, or both. "
            "Void subtraction is applied before splitting."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    radius = args.diameter / 2.0
    base_body_meshes = build_sphere_body_meshes(
        radius_mm=radius,
        lat_segments=args.lat_segments,
        lon_segments=args.lon_segments,
        void_shape=args.void_shape,
        void_size_mm=args.void_size,
        split_half=args.split_half,
    )
    sphere_meshes = tile_body_meshes_2d(
        base_meshes=base_body_meshes,
        array_x=args.array_x,
        array_y=args.array_y,
        array_spacing_mm=args.array_spacing,
    )
    anchor_meshes = build_corner_anchor_meshes(
        enabled=args.corner_anchors,
        sphere_diameter_mm=args.diameter,
        sphere_radius_mm=radius,
        array_x=args.array_x,
        array_y=args.array_y,
        array_spacing_mm=args.array_spacing,
    )
    all_meshes = sphere_meshes + anchor_meshes
    triangles = flatten_meshes(all_meshes)
    output_format = resolve_output_format(args.format, args.output)

    if args.split_half == "up":
        solid_name = "sphere_upper_half"
    elif args.split_half == "down":
        solid_name = "sphere_lower_half"
    elif args.split_half == "both":
        solid_name = "sphere_split_halves"
    elif args.void_shape != "none":
        solid_name = "sphere_with_internal_void"
    else:
        solid_name = "sphere_solid"

    mesh_list: List[List[Triangle]] | None = None
    if output_format == "3mf":
        mesh_list = list(all_meshes)

    write_mesh_output(
        args.output, output_format, solid_name, triangles, mesh_list=mesh_list
    )

    if args.void_shape == "none":
        void_desc = "internal_void=none"
    else:
        void_desc = f"internal_void={args.void_shape}({args.void_size} mm)"
    split_desc = f"split_half={args.split_half}"
    anchor_desc = "corner_anchors=on" if args.corner_anchors else "corner_anchors=off"
    if args.array_x == 1 and args.array_y == 1:
        array_desc = "array=1x1"
    else:
        array_desc = (
            f"array={args.array_x}x{args.array_y}, "
            f"spacing={args.array_spacing} mm"
        )
    sphere_instances = args.array_x * args.array_y
    sphere_bodies = len(sphere_meshes)
    anchor_bodies = len(anchor_meshes)
    total_bodies = len(all_meshes)
    print(
        f"Generated '{args.output}' with diameter={args.diameter} mm, "
        f"format={output_format}, {split_desc}, {void_desc}, {array_desc}, "
        f"{anchor_desc}, sphere_instances={sphere_instances}, "
        f"sphere_bodies={sphere_bodies}, anchor_bodies={anchor_bodies}, "
        f"bodies={total_bodies}, facets={len(triangles)}."
    )


if __name__ == "__main__":
    main()
