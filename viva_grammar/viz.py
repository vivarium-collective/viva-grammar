"""Interactive visualization of the hypoxia -> motility reference study.

Renders a self-contained HTML page (Plotly, CDN-loaded) with:
  1. an **animated spatial view** — the oxygen field (single-hue blue ramp) with
     the cell cohort overlaid, coloured by state (tumor / motile tumor / necrotic),
     a play/pause control and a time slider; and
  2. a **dynamics dashboard** — motile-cell fraction over time with the hypoxic and
     reoxygenation phases shaded, the population stack, and the radial-distribution
     traces that evidence the invasive plume and the central necrotic core.

Palette: oxygen = sequential Blues; cell states categorical slate / orange /
near-black (legend + marker size as secondary encoding).  Dark-mode aware via a
prefers-color-scheme block around the Plotly template.
"""
from __future__ import annotations

import json
from typing import Dict, List

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# categorical cell-state colors (see module docstring)
C_TUMOR = "#6b7280"
C_MOTILE = "#e8663c"
C_NECROTIC = "#111827"
PLOTLY_CDN = "https://cdn.plot.ly/plotly-3.1.0.min.js"


def _frac(pops: Dict[str, float]) -> float:
    live = pops.get("live_tumor", 0.0) + pops.get("live_motile_tumor", 0.0)
    return (pops.get("live_motile_tumor", 0.0) / live) if live else 0.0


def _round_grid(grid):
    """Round the oxygen grid to 1 decimal to keep the committed HTML small."""
    return [[round(float(v), 1) for v in row] for row in grid]


def _spatial_figure(frames: List[dict], domain: dict) -> go.Figure:
    x0, x1 = domain["x_min"], domain["x_max"]
    y0, y1 = domain["y_min"], domain["y_max"]

    def cell_traces(fr):
        groups = {"tumor": ([], []), "motile tumor": ([], []), "necrotic": ([], [])}
        for (x, y), st in zip(fr["cell_xy"], fr["cell_state"]):
            key = "necrotic" if st == "necrotic" else st
            if key not in groups:
                key = "tumor"
            groups[key][0].append(x); groups[key][1].append(y)
        specs = [
            ("tumor", C_TUMOR, 7), ("motile tumor", C_MOTILE, 8), ("necrotic", C_NECROTIC, 5),
        ]
        traces = []
        for name, color, size in specs:
            xs, ys = groups[name]
            traces.append(go.Scatter(
                x=xs, y=ys, mode="markers", name=name,
                marker=dict(color=color, size=size, line=dict(width=0.5, color="rgba(255,255,255,0.6)"),
                            opacity=0.55 if name == "necrotic" else 0.95),
                hovertemplate=f"{name}<br>x=%{{x:.0f}} y=%{{y:.0f}}<extra></extra>",
            ))
        return traces

    f0 = frames[0]
    heat0 = go.Heatmap(
        z=_round_grid(f0["oxygen"]), x0=x0, dx=(x1 - x0) / len(f0["oxygen"][0]),
        y0=y0, dy=(y1 - y0) / len(f0["oxygen"]),
        colorscale="Blues", zmin=0, zmax=38, opacity=0.85,
        colorbar=dict(title="O₂ (mmHg)", len=0.75, thickness=12),
        hovertemplate="O₂=%{z:.1f} mmHg<extra></extra>",
    )
    fig = go.Figure(data=[heat0, *cell_traces(f0)])

    plotly_frames = []
    for fr in frames:
        heat = go.Heatmap(z=_round_grid(fr["oxygen"]), x0=x0, dx=(x1 - x0) / len(fr["oxygen"][0]),
                          y0=y0, dy=(y1 - y0) / len(fr["oxygen"]),
                          colorscale="Blues", zmin=0, zmax=38, opacity=0.85, showscale=True,
                          colorbar=dict(title="O₂ (mmHg)", len=0.75, thickness=12))
        plotly_frames.append(go.Frame(data=[heat, *cell_traces(fr)], name=f"{fr['time']:.0f}"))
    fig.frames = plotly_frames

    steps = [dict(method="animate", label=f"{fr['time']:.0f}",
                  args=[[f"{fr['time']:.0f}"],
                        dict(mode="immediate", frame=dict(duration=0, redraw=True),
                             transition=dict(duration=0))]) for fr in frames]
    fig.update_layout(
        title="Oxygen field & cell cohort over time (min)",
        xaxis=dict(title="x (µm)", range=[x0, x1], constrain="domain"),
        yaxis=dict(title="y (µm)", range=[y0, y1], scaleanchor="x", scaleratio=1),
        legend=dict(orientation="h", y=1.02, x=0, bgcolor="rgba(0,0,0,0)"),
        margin=dict(l=60, r=20, t=90, b=50),
        updatemenus=[dict(type="buttons", x=0, y=1.14, xanchor="left", showactive=False,
                          buttons=[
                              dict(label="▶ Play", method="animate",
                                   args=[None, dict(frame=dict(duration=350, redraw=True),
                                                    fromcurrent=True, transition=dict(duration=0))]),
                              dict(label="⏸ Pause", method="animate",
                                   args=[[None], dict(mode="immediate",
                                                      frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, y=0, x=0.08, len=0.9, pad=dict(t=40),
                      currentvalue=dict(prefix="t = ", suffix=" min"), steps=steps)],
    )
    return fig


def _dynamics_figure(frames: List[dict], reoxy_time: float) -> go.Figure:
    t = [fr["time"] for fr in frames]
    tumor = [fr["populations"].get("live_tumor", 0) for fr in frames]
    motile = [fr["populations"].get("live_motile_tumor", 0) for fr in frames]
    necro = [fr["populations"].get("necrotic", 0) for fr in frames]
    frac = [_frac(fr["populations"]) for fr in frames]
    mot_r = [fr["radial"].get("mean_motile_radius", 0) for fr in frames]
    nonmot_r = [fr["radial"].get("mean_nonmotile_radius", 0) for fr in frames]
    necro_r = [fr["radial"].get("mean_necrotic_radius", 0) for fr in frames]

    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07,
        subplot_titles=("Motile-cell fraction (phenotype induction & reversion)",
                        "Population by state", "Mean radial distance from center (plume & core)"),
    )
    # phase shading: reoxygenation onward
    tmax = t[-1] if t else 1
    fig.add_vrect(x0=reoxy_time, x1=tmax, fillcolor="rgba(78,121,167,0.10)", line_width=0,
                  annotation_text="reoxygenation", annotation_position="top left", row="all", col=1)

    fig.add_trace(go.Scatter(x=t, y=frac, name="motile fraction", mode="lines",
                             line=dict(color=C_MOTILE, width=2),
                             hovertemplate="t=%{x:.0f}<br>motile frac=%{y:.3f}<extra></extra>"), row=1, col=1)
    fig.add_hline(y=0.15, line=dict(color="rgba(0,0,0,0.3)", dash="dot"), row=1, col=1)

    fig.add_trace(go.Scatter(x=t, y=tumor, name="tumor", mode="lines", stackgroup="pop",
                             line=dict(width=0.5, color=C_TUMOR), fillcolor="rgba(107,114,128,0.6)"), row=2, col=1)
    fig.add_trace(go.Scatter(x=t, y=motile, name="motile tumor", mode="lines", stackgroup="pop",
                             line=dict(width=0.5, color=C_MOTILE), fillcolor="rgba(232,102,60,0.7)"), row=2, col=1)
    fig.add_trace(go.Scatter(x=t, y=necro, name="necrotic", mode="lines", stackgroup="pop",
                             line=dict(width=0.5, color=C_NECROTIC), fillcolor="rgba(17,24,39,0.6)"), row=2, col=1)

    fig.add_trace(go.Scatter(x=t, y=mot_r, name="motile radius", mode="lines",
                             line=dict(color=C_MOTILE, width=2)), row=3, col=1)
    fig.add_trace(go.Scatter(x=t, y=nonmot_r, name="non-motile radius", mode="lines",
                             line=dict(color=C_TUMOR, width=2)), row=3, col=1)
    fig.add_trace(go.Scatter(x=t, y=necro_r, name="necrotic radius", mode="lines",
                             line=dict(color=C_NECROTIC, width=2, dash="dash")), row=3, col=1)

    fig.update_yaxes(title_text="fraction", row=1, col=1, rangemode="tozero")
    fig.update_yaxes(title_text="cells", row=2, col=1, rangemode="tozero")
    fig.update_yaxes(title_text="µm", row=3, col=1, rangemode="tozero")
    fig.update_xaxes(title_text="time (min)", row=3, col=1)
    fig.update_layout(margin=dict(l=60, r=20, t=50, b=50),
                      legend=dict(orientation="h", y=-0.12), hovermode="x unified")
    return fig


_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Hypoxia → Motility (grammar)</title>
<script src="{cdn}"></script>
<style>
  :root {{ --bg:#fcfcfb; --ink:#1a1a1a; --muted:#5b5b5b; --card:#ffffff; --border:#e7e7e3; }}
  :root:not([data-theme="light"]) {{ }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{ --bg:#15171a; --ink:#f0f0ee; --muted:#a7a7a2; --card:#1d2024; --border:#2c2f34; }}
  }}
  :root[data-theme="dark"] {{ --bg:#15171a; --ink:#f0f0ee; --muted:#a7a7a2; --card:#1d2024; --border:#2c2f34; }}
  html,body {{ background:var(--bg); color:var(--ink); margin:0;
    font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
  .wrap {{ max-width:1100px; margin:0 auto; padding:24px 16px 64px; }}
  h1 {{ font-size:1.5rem; margin:0 0 4px; }}
  p.sub {{ color:var(--muted); margin:0 0 20px; max-width:70ch; }}
  .card {{ background:var(--card); border:1px solid var(--border); border-radius:12px;
    padding:12px; margin:18px 0; }}
  .bands {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr)); gap:10px; margin:18px 0; }}
  .band {{ background:var(--card); border:1px solid var(--border); border-radius:10px; padding:12px 14px; }}
  .band .v {{ font-size:1.25rem; font-weight:650; }}
  .band .ok {{ color:#2e7d46; }} .band .no {{ color:#b23b3b; }}
  .band .lbl {{ color:var(--muted); font-size:.85rem; }}
  .plot {{ width:100%; }}
</style></head>
<body><div class="wrap">
  <h1>Hypoxia → motility, encoded as a cell-behavior grammar</h1>
  <p class="sub">The verbatim PhysiCell hypoxia rules, executed by the generic
  <code>GrammarRulesProcess</code> on a minimal oxygen-field + point-agent reference
  substrate. Low oxygen drives a transition to a motile phenotype that disseminates
  from the hypoxic core and reverts on reoxygenation — the paper's transient
  post-hypoxic migration.</p>
  <div class="bands">{bands}</div>
  <div class="card"><div id="spatial" class="plot"></div></div>
  <div class="card"><div id="dynamics" class="plot"></div></div>
</div>
<script>
  const spatial = {spatial};
  const dynamics = {dynamics};
  const cfg = {{responsive:true, displayModeBar:false}};
  // Pass the whole figure object (data + layout + frames) so the slider/play work.
  Plotly.newPlot('spatial', spatial, cfg);
  Plotly.newPlot('dynamics', dynamics, cfg);
</script>
</body></html>
"""


def _band_card(label, value, ok, detail):
    cls = "ok" if ok else "no"
    mark = "✓" if ok else "✗"
    return (f'<div class="band"><div class="lbl">{label}</div>'
            f'<div class="v {cls}">{mark} {value}</div>'
            f'<div class="lbl">{detail}</div></div>')


def render_hypoxia_html(frames: List[dict], domain: dict, reoxy_time: float, bands: dict) -> str:
    spatial = _spatial_figure(frames, domain)
    dynamics = _dynamics_figure(frames, reoxy_time)
    band_html = "".join(
        _band_card(b["label"], b["value"], b["ok"], b["detail"]) for b in bands["cards"]
    )
    # Build frames JSON for the spatial animation (Plotly.newPlot takes data/layout; frames set after).
    spatial_json = spatial.to_plotly_json()
    import plotly.io as pio
    spatial_dict = json.loads(pio.to_json(spatial))
    dynamics_dict = json.loads(pio.to_json(dynamics))
    return _PAGE.format(
        cdn=PLOTLY_CDN,
        bands=band_html,
        spatial=json.dumps(spatial_dict),
        dynamics=json.dumps(dynamics_dict),
    )
