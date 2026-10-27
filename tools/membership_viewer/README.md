# Reusable membership inspector

The inspector reads the saved event memberships produced by the joint-surface
campaign. It never fits a surface or changes a selection. Both views show measured
focal-plane position, measured focal-plane slope, and event delta. Positions are
in cm, slopes are displayed in mrad, and delta is in percent.

## Native window on ifarm

Connect with X forwarding (`ssh -Y`) and a running X server on your local computer.
From the repository root, after the campaign label stage:

```sh
python3 -m joint_surface.viewer /path/to/CAMPAIGN
```

This opens a Matplotlib/Tk window; it uses no browser or web server. It requires
NumPy, Matplotlib, and Python's Tk interface. Tk is usually a system package
(`python3-tk` on Debian/Ubuntu); installing Matplotlib with pip alone does not
provide it. If the window cannot open, check `echo "$DISPLAY"` and X forwarding.
With a slow forwarded display, reduce the plotting budget:

```sh
python3 -m joint_surface.viewer /path/to/CAMPAIGN --max-points 4000
```

Select the foil/run, saved configuration, delta interval, and X/Y color rule from
the menus. Color either plane by its own membership, or color both planes by the
same X row or Y column. Crosses show rejected events, with separate legend entries
for each saved rejection reason. The event filter can isolate events accepted in
both axes or rejected in either axis. Drag to rotate. Use the toolbar to pan,
zoom, and save an image. Click a point to highlight the same entry in both planes
and read its X/Y membership, rejection reasons, and saved training-mask metadata.
Click a legend item to hide/show that membership; shared X or Y groups are linked.

The native window uses the training mask only as information about a saved event.
It does not train or count events during viewing. Optional comparison labels are
also display information and are not required.

## Self-contained offline Firefox alternative

Install the optional Plotly Python dependency once in the analysis environment:

```sh
python3 -m pip install plotly
python3 -m joint_surface.viewer /path/to/CAMPAIGN \
  --html /path/to/CAMPAIGN/membership.html
firefox /path/to/CAMPAIGN/membership.html
```

The export embeds the installed Plotly JavaScript and all plotted event records in
one HTML file. Opening that file needs neither internet access nor a server. It
can also be copied to another machine. Alternatively supply an existing local
bundle with `--plotly-js /path/to/plotly.min.js`. The program never downloads a
JavaScript library. Use `--force` to intentionally replace a previous HTML file.
Firefox must permit WebGL for Plotly's 3D plots. If Firefox cannot create a WebGL
context over X forwarding, use the native window or copy the HTML to a machine
with local graphics.

Hover for event IDs and memberships. Click a point to highlight it in both planes.
Click a legend item to hide a group, or double-click to isolate it. Shared X/Y
colors link legend visibility between planes. No data are uploaded.

## Population counts and sampling

By default, the display contains at most 12,000 events per foil. Selection is
deterministic: up to one quarter of this budget is shared between rejection
categories, and remaining places are sampled randomly from the remaining events.
The same displayed events are used across configurations. This deliberately
emphasizes some rare rejections; the plot is not an unbiased density estimate.

The full interval counts and plotted counts are displayed separately. The native
status line reports the sample size before the optional accepted/rejected event
filter; the HTML reports the number after that filter. Events with a nonfinite
coordinate cannot be plotted but remain in the full population counts. The
delta intervals include their lower bound and exclude their upper bound. They
filter saved events for inspection and never refit a surface. `--max-points 0`
plots every event with finite coordinates; use this cautiously for large replays.

## Saved-data contract

The default manifest is `CAMPAIGN/00_joint_surface_work/viewer_manifest.json`.
An explicit manifest or its containing directory can also be passed to the CLI.
Schema version 1 contains a `cases` list. Every case has a unique `key`, a `label`,
an `events` NPZ path, a mapping `predictions` of saved configuration names to NPZ
paths, and optional `delta_edges` (default `[-10,-8,-5,0,5,10]`). Paths are relative
to the manifest directory; an absolute path is accepted. A single `npz` path can
stand in for both `events` and the sole `main` prediction.

Event NPZ arrays are `entry` (exact int64 original replay entry), `xfp`, `xpfp`,
`yfp`, `ypfp`, and `delta`. All arrays are one-dimensional in identical event
order. Stored positions are cm, stored slopes radians, and delta percent.
Prediction arrays are `x_label`, `y_label`, `accepted_x`, `accepted_y`, `reason_x`,
and `reason_y`. Rejected indices may be -1. An accepted mask is authoritative:
a proposed nonnegative index does not by itself make a rejected event accepted.
If masks are absent, nonnegative labels mean accepted. `physical_x/y` or
`latent_x/y` are accepted as label aliases. Optional one-dimensional arrays with
prefixes `proposed_`, `latent_`, `train_`, `reference_`, or `final_` appear in event
details. Arrays of Python objects/pickled data are not accepted.

The prediction `entry` array, when supplied, must exactly match event order.
No event join, sorting, label reconciliation, or reconstruction is performed by
the viewer. Different saved configurations can be added to the manifest without
changing or rebuilding the interface.

Targeted checks:

```sh
MPLBACKEND=Agg python3 -m unittest discover -s tests -p test_joint_surface_viewer.py -v
```
