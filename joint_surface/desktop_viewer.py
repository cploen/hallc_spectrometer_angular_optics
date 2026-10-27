"""Native Tk/Matplotlib inspector. No browser, web server, or network is used."""
from __future__ import annotations

import numpy as np

from .viewer import counts, event_text, trace_groups


def picked_index(case, indices, axes, mouse_event):
    """Find the nearest event in screen coordinates after a 3D rotation.

    Matplotlib sorts a 3D scatter by depth while drawing; its pick indices can
    therefore differ from the input event order. Project the original measured
    coordinates instead of treating a depth-sorted pick index as an event ID.
    """
    from mpl_toolkits.mplot3d import proj3d
    plane = 0 if axes.get_xlabel().startswith("xfp") else 2
    xyz = case.coordinates[indices]
    x, y, _ = proj3d.proj_transform(xyz[:, plane], xyz[:, plane + 1], xyz[:, 4], axes.get_proj())
    pixels = axes.transData.transform(np.column_stack((x, y)))
    distance2 = (pixels[:, 0] - mouse_event.x)**2 + (pixels[:, 1] - mouse_event.y)**2
    return int(indices[np.argmin(distance2)])


def draw_case(figure, case, config, interval=None, color_mode="each", event_filter="all", cameras=None):
    """Render measured coordinates; also usable with Agg for a saved screenshot."""
    figure.clear()
    axes, artists, legends = [], {}, {}
    for number, axis in enumerate(("x", "y")):
        ax = figure.add_subplot(1, 2, number + 1, projection="3d")
        axes.append(ax)
        j = number * 2
        for name, color, marker, indices in trace_groups(case, config, axis, color_mode, interval, event_filter):
            xyz = case.coordinates[indices]
            artist = ax.scatter(xyz[:, j], xyz[:, j + 1], xyz[:, 4], s=7 if marker == "o" else 18,
                                c=color, marker=marker, alpha=.68 if marker == "o" else .95,
                                label=f"{name} ({len(indices):,} shown)", depthshade=False, picker=5)
            artists[artist] = (number, name, indices)
        ax.set_xlabel(axis + "fp (cm)")
        ax.set_ylabel(axis + "pfp (mrad)")
        ax.set_zlabel("delta (%)")
        ax.set_title(axis.upper() + " focal plane")
        ax.view_init(elev=15, azim=-60)
        if cameras:
            ax.view_init(elev=cameras[number][0], azim=cameras[number][1])
        ax.set_box_aspect((1, 1, 1.1))
        if artists:
            legend = ax.legend(loc="upper left", fontsize=7, markerscale=2, handletextpad=.3,
                               borderpad=.3, framealpha=.85)
            handles = getattr(legend, "legend_handles", getattr(legend, "legendHandles", []))
            local = [(a, info) for a, info in artists.items() if info[0] == number]
            for handle, text, (artist, info) in zip(handles, legend.get_texts(), local):
                for target in (handle, text):
                    target.set_picker(True)
                    legends[target] = info[1]
    figure.subplots_adjust(left=.02, right=.98, bottom=.06, top=.95, wspace=.03)
    return axes, artists, legends


def launch(cases):
    try:
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure
    except ImportError as exc:
        raise ValueError("The desktop viewer needs Matplotlib and Python Tk support. On Linux install the distribution's python3-tk package, or use --html with Firefox.") from exc
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise ValueError("Cannot open the native display. Connect to ifarm with ssh -Y and a running X server; check DISPLAY. A headless shell can export --html for offline Firefox inspection.") from exc
    root.title("Joint surface membership — measured focal-plane coordinates")
    root.geometry("1450x950")
    root.minsize(1000, 650)
    controls = ttk.Frame(root, padding=(8, 8))
    controls.pack(fill="x")
    selection = {"case": 0, "config": next(iter(cases[0].configurations)), "interval": None,
                 "color": "each", "filter": "all"}
    variables = {name: tk.StringVar() for name in ("case", "config", "interval", "color", "filter")}
    boxes = {}
    labels = {"case": "Foil / run", "config": "Saved configuration", "interval": "Delta interval",
              "color": "Color membership", "filter": "Events"}
    for col, name in enumerate(labels):
        frame = ttk.Frame(controls)
        frame.grid(row=0, column=col, sticky="ew", padx=(0, 10))
        ttk.Label(frame, text=labels[name]).pack(anchor="w")
        box = ttk.Combobox(frame, state="readonly", textvariable=variables[name], width=35 if name == "case" else 23)
        box.pack(fill="x")
        boxes[name] = box
        controls.columnconfigure(col, weight=2 if name == "case" else 1)
    boxes["case"]["values"] = [c.label for c in cases]
    boxes["color"]["values"] = ["Each plane's membership", "X row in both planes", "Y column in both planes"]
    boxes["filter"]["values"] = ["All events", "Both axes accepted", "Either axis rejected"]
    boxes["case"].current(0)
    boxes["color"].current(0)
    boxes["filter"].current(0)
    status = tk.StringVar()
    ttk.Label(root, textvariable=status, padding=(8, 0), wraplength=1400).pack(fill="x")
    ttk.Label(root, text="Drag a plot to rotate. Use the toolbar to pan, zoom, or save. Click an event for its entry ID; click a legend item to hide/show that membership.", padding=(8, 3)).pack(fill="x")
    figure = Figure(figsize=(14, 8), dpi=100)
    canvas = FigureCanvasTkAgg(figure, root)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    NavigationToolbar2Tk(canvas, root)
    detail = tk.StringVar(value="Click a point to inspect the same event in both focal planes.")
    ttk.Label(root, textvariable=detail, wraplength=1400, padding=8).pack(fill="x")
    plot_state = {"axes": [], "artists": {}, "legends": {}, "highlights": []}

    def update_options():
        c = cases[selection["case"]]
        configs = list(c.configurations)
        boxes["config"]["values"] = configs
        if selection["config"] not in configs:
            selection["config"] = configs[0]
        boxes["config"].current(configs.index(selection["config"]))
        boxes["interval"]["values"] = ["All delta", *[f"{lo:g}% ≤ delta < {hi:g}%" for lo, hi in zip(c.delta_edges[:-1], c.delta_edges[1:])]]
        if selection["interval"] is not None and selection["interval"] >= len(c.delta_edges)-1:
            selection["interval"] = None
        boxes["interval"].current(0 if selection["interval"] is None else selection["interval"] + 1)

    def redraw(reset=False):
        c = cases[selection["case"]]
        cameras = None if reset else [(a.elev, a.azim) for a in plot_state["axes"]] or None
        args = (figure, c, selection["config"], selection["interval"], selection["color"], selection["filter"], cameras)
        axes, artists, legends = draw_case(*args)
        plot_state.update(axes=axes, artists=artists, legends=legends, highlights=[])
        full = counts(c, selection["config"], selection["interval"])
        shown = counts(c, selection["config"], selection["interval"], shown=True)
        status.set(f"Full interval population: {full['events']:,}; X accepted {full['x']:,}; Y accepted {full['y']:,}; both {full['both']:,}. "
                   f"Sample before event filter: {shown['events']:,} plotted from {len(c.entry):,} full-case events. "
                   "Display sampling preserves some rejection categories; plotted density is not a population density estimate.")
        canvas.draw_idle()

    def change(name):
        if name == "case":
            selection["case"] = boxes[name].current()
            update_options()
        elif name == "config":
            selection["config"] = variables[name].get()
        elif name == "interval":
            selection["interval"] = None if boxes[name].current() == 0 else boxes[name].current() - 1
        elif name == "color":
            selection["color"] = ("each", "x", "y")[boxes[name].current()]
        else:
            selection["filter"] = ("all", "accepted", "rejected")[boxes[name].current()]
        redraw()

    for name, box in boxes.items():
        box.bind("<<ComboboxSelected>>", lambda event, name=name: change(name))

    def pick(event):
        if event.artist in plot_state["legends"]:
            name = plot_state["legends"][event.artist]
            matching = [a for a, info in plot_state["artists"].items() if info[1] == name]
            visible = not matching[0].get_visible()
            for a in matching:
                a.set_visible(visible)
            for a, label in plot_state["legends"].items():
                if label == name:
                    a.set_alpha(1 if visible else .25)
            canvas.draw_idle()
            return
        if event.artist not in plot_state["artists"] or not len(event.ind):
            return
        plane, _, indices = plot_state["artists"][event.artist]
        c = cases[selection["case"]]
        index = picked_index(c, indices, plot_state["axes"][plane], event.mouseevent)
        detail.set(event_text(c, selection["config"], index))
        for artist in plot_state["highlights"]:
            artist.remove()
        plot_state["highlights"] = []
        for j, ax in enumerate(plot_state["axes"]):
            xyz = c.coordinates[index]
            h = ax.scatter([xyz[2*j]], [xyz[2*j+1]], [xyz[4]], s=85,
                           facecolors="none", edgecolors="black", linewidths=1.5, depthshade=False)
            plot_state["highlights"].append(h)
        canvas.draw_idle()

    canvas.mpl_connect("pick_event", pick)
    update_options()
    redraw(reset=True)
    root.mainloop()
