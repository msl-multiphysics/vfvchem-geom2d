"""PyVista rendering of domains, boundaries, and interfaces."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


def _label_colors(regions: list, cmap_name: str) -> dict:
    """RGB color per label, in order of first appearance."""
    cmap = plt.get_cmap(cmap_name)
    labels = dict.fromkeys(label for label, _ in regions)
    return {label: cmap(i % cmap.N)[:3] for i, label in enumerate(labels)}


def _contrast_text_color(rgb) -> str:
    r, g, b = rgb[:3]
    # Rec. 601 luma; light backgrounds get black text
    return "black" if (0.299 * r + 0.587 * g + 0.114 * b) > 0.55 else "white"


def _domain_surface(loops: list, z: float):
    """Triangulated domain filling the outer loops and leaving the holes empty."""
    from vtkmodules.vtkFiltersGeneral import vtkContourTriangulator

    pts = []
    lines = []
    for lp in loops:
        # loops repeat their first point
        n = len(pts)
        pts.extend((x, y, z) for x, y in lp[:-1])
        lines.extend([len(lp), *range(n, len(pts)), n])
    tri = vtkContourTriangulator()
    tri.SetInputData(pv.PolyData(np.asarray(pts, dtype=float), lines=lines))
    tri.Update()
    return pv.wrap(tri.GetOutput())


def _text_button_texture(label: str, font_path: str | None = None, font_size: int = 16):
    """Build a PyVista ImageData texture with a labeled button appearance."""
    from PIL import Image, ImageDraw, ImageFont

    if font_path:
        font = ImageFont.truetype(font_path, font_size)
    else:
        try:
            font = ImageFont.truetype("Arial.ttf", font_size)
        except OSError:
            font = ImageFont.load_default()

    probe = Image.new("RGB", (1, 1))
    draw = ImageDraw.Draw(probe)
    bbox = draw.textbbox((0, 0), label, font=font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    pad_x, pad_y = 14, 8
    width = text_w + 2 * pad_x
    height = text_h + 2 * pad_y

    img = Image.new("RGB", (width, height), (235, 235, 235))
    draw = ImageDraw.Draw(img)
    draw.rectangle([0, 0, width - 1, height - 1], outline=(60, 60, 60), width=2)
    draw.text((pad_x - bbox[0], pad_y - bbox[1]), label, fill=(20, 20, 20), font=font)

    arr = np.asarray(img, dtype=np.uint8)  # (height, width, 3), y=0 at top
    # VTK ImageData is x-fastest with y=0 at bottom
    arr = np.flipud(arr)
    button = pv.ImageData(dimensions=(width, height, 1))
    button.point_data["texture"] = arr.reshape(width * height, 3)
    return button, width, height


def _add_reset_camera_button(plotter, *, margin: float = 10.0):
    """Add a clickable 'Reset Camera' text button at the bottom-right."""
    from vtkmodules.vtkCommonCore import vtkCommand
    from vtkmodules.vtkInteractionWidgets import (
        vtkButtonWidget,
        vtkTexturedButtonRepresentation2D,
    )

    texture, width, height = _text_button_texture(
        "Reset Camera",
        font_path=_resolve_label_font("IBM Plex Sans"),
    )
    win_w, win_h = plotter.window_size
    position = (win_w - width - margin, margin)

    # Same texture for both states so the control does not look like a toggle
    button_rep = vtkTexturedButtonRepresentation2D()
    button_rep.SetNumberOfStates(2)
    button_rep.SetButtonTexture(0, texture)
    button_rep.SetButtonTexture(1, texture)
    button_rep.SetState(0)
    button_rep.SetPlaceFactor(1)
    button_rep.PlaceWidget(
        [position[0], position[0] + width, position[1], position[1] + height, 0.0, 0.0]
    )

    button_widget = vtkButtonWidget()
    button_widget.SetInteractor(plotter.iren.interactor)
    button_widget.SetRepresentation(button_rep)
    button_widget.SetCurrentRenderer(plotter.renderer)
    button_widget.On()

    def _on_click(_widget, _event):
        plotter.view_xy()
        plotter.reset_camera()
        button_rep.SetState(0)

    button_widget.AddObserver(vtkCommand.StateChangedEvent, _on_click)
    plotter.widgets.button_widgets.append(button_widget)
    return button_widget


def _resolve_label_font(family: str = "IBM Plex Sans") -> str | None:
    """Return a TTF path for the label font, or None if unavailable."""
    from matplotlib import font_manager

    try:
        path = font_manager.findfont(
            font_manager.FontProperties(family=family),
            fallback_to_default=False,
        )
    except (ValueError, RuntimeError):
        return None
    if path and "DejaVu" not in path:  # matplotlib's silent fallback name
        return path
    # Explicit scan in case findfont fell back quietly on some versions
    for entry in font_manager.fontManager.ttflist:
        if entry.name == family and entry.style == "normal" and int(entry.weight) == 400:
            return entry.fname
    for entry in font_manager.fontManager.ttflist:
        if entry.name == family and entry.style == "normal":
            return entry.fname
    return None


def render_geometry(
    title: str,
    z: float,
    domains: list,
    boundaries: list,
    interfaces: list,
    *,
    show_domains: bool = True,
    show_boundaries: bool = True,
    show_interfaces: bool = True,
    show_labels: bool = True,
):
    """Open an interactive PyVista window with the geometry.

    domains are (label, [closed loop, ...]) with ccw outer loops and cw holes.
    boundaries and interfaces are (label, [polyline, ...]). Points are (x, y).
    Regions with the same label share a color.
    """
    plotter = pv.Plotter(title=title)
    plotter.set_background("white")
    plotter.show_axes()

    label_items = []  # (label points, label, color)

    # domain fill only; the triangulation is not a mesh, so its edges are hidden
    if show_domains and domains:
        colors = _label_colors(domains, "Pastel1")
        for label, loops in domains:
            surf = _domain_surface(loops, z)
            plotter.add_mesh(surf, color=colors[label], lighting=False)
            areas = surf.compute_cell_sizes(length=False, volume=False)["Area"]
            label_items.append((surf.cell_centers().points[[np.argmax(areas)]], label, colors[label]))

    # vertex markers show how curves are approximated by segments
    for regions, cmap_name, show in (
        (boundaries, "Set1", show_boundaries),
        (interfaces, "Set2", show_interfaces),
    ):
        if not show or not regions:
            continue
        colors = _label_colors(regions, cmap_name)
        for label, lines in regions:
            for line in lines:
                pts = np.asarray([(x, y, z) for x, y in line], dtype=float)
                plotter.add_mesh(pv.lines_from_points(pts), color=colors[label], line_width=4)
                plotter.add_points(pts, color=colors[label], point_size=8, render_points_as_spheres=True)
                i = (len(pts) - 1) // 2
                label_items.append((0.5 * (pts[[i]] + pts[[i + 1]]), label, colors[label]))

    if show_labels:
        font_file = _resolve_label_font("IBM Plex Sans")
        for pts, label, rgb in label_items:
            plotter.add_point_labels(
                pts,
                [label] * len(pts),
                point_size=5,
                font_size=14,
                bold=False,
                font_file=font_file,
                text_color=_contrast_text_color(rgb),
                shape_color=rgb,
                shape_opacity=0.8,
                always_visible=True,
            )

    plotter.view_xy()

    _add_reset_camera_button(plotter)

    plotter.show()
