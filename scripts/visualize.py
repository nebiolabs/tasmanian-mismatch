#!/usr/bin/env python3
"""
Bokeh interactive visualization for tasmanian-mismatch output.

Usage:
    python scripts/visualize.py output.tsv                  # opens in browser
    python scripts/visualize.py output.tsv -o plot.html      # saves standalone HTML

Reads the 5-column TSV produced by tasmanian-mismatch:
    base_change  read_num  reference_order  read_position|fragment_position  count

Features:
  - One colored line per mismatch type (A>C, A>G, ... excludes self-matches)
  - Hover tooltip showing position, count, base_change
  - Read1 / Read2 toggle widgets (hidden when input is Insert mode)
    - Optional normalization: X>Y / (X>A + X>C + X>G + X>T)
  - Plot occupies 3/4 of the viewport height
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Literal

import pandas as pd
from bokeh.embed import file_html
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    CheckboxGroup,
    ColumnDataSource,
    CustomJS,
    Div,
    HoverTool,
    Legend,
    LegendItem,
    Spacer,
    InlineStyleSheet,
)
from bokeh.palettes import Category20_16
from bokeh.plotting import figure
from bokeh.resources import Resources


MISMATCHES = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]

PALETTE = {m: Category20_16[i] for i, m in enumerate(MISMATCHES)}

READ_DASH_PATTERNS = {1: "solid", 2: "dashed"}
READ_MARKER_SHAPES = {1: "circle", 2: "triangle"}
CHECKBOX_TEXT_STYLE: dict[str, str | None] = {"font-size": "16pt", "line-height": "1.6"}
PANEL_TITLE_STYLE: dict[str, str | None] = {"font-size": "18pt", "font-weight": "bold", "margin-bottom": "8px"}

PAGE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{{ title }}</title>
  {{ bokeh_css }}
  {{ bokeh_js }}
  <style>
    html, body { margin: 0; padding: 0; width: 100%%; height: 100%%; }
    #plot-container { width: 95%%; height: 75vh; margin: 1em auto; }
  </style>
</head>
<body>
  <div id="plot-container">
    {{ plot_div }}
  </div>
  {{ plot_script }}
</body>
</html>
"""


def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    pos_col = [c for c in df.columns if c in ("read_position", "fragment_position")]
    if not pos_col:
        sys.exit(f"Expected 'read_position' or 'fragment_position' column. Got: {list(df.columns)}")
    df = df.rename(columns={pos_col[0]: "position"})
    df["is_insert_mode"] = pos_col[0] == "fragment_position"
    return df


def empty_plot_data() -> dict:
    return {
        "position": [],
        "y": [],
        "y_count": [],
        "y_norm": [],
        "base_change": [],
        "read_num": [],
        "reference_order": [],
    }


def make_renderer_pair(
    p,
    source,
    color: str,
    dash: str,
    marker: str,
    line_alpha: float,
    point_alpha: float,
    visible: bool,
):
    line = p.line(
        "position", "y", source=source,
        line_width=8, color=color, alpha=line_alpha,
        line_dash=dash,
        visible=visible,
    )
    points = p.scatter(
        "position", "y", source=source,
        size=4, color=color, alpha=point_alpha,
        marker=marker,
        visible=visible,
    )
    bars = p.vbar(
        x="position",
        top="y",
        width=0.8,
        source=source,
        fill_color=color,
        line_color=color,
        fill_alpha=point_alpha,
        line_alpha=line_alpha,
        visible=False,
    )
    return line, points, bars


def build_button_text_css(font_px: int = 28) -> str:
    return (
        f":host {{ font-size: {font_px}px !important; font-weight: 700 !important; }}"
        f":host button, :host .bk-btn {{ font-size: {font_px}px !important; font-weight: 700 !important; border-radius: 20px !important; line-height: 1.1 !important; }}"
        f":host span, :host .bk-btn-label, :host .bk-btn-text {{ font-size: {font_px}px !important; font-weight: 700 !important; line-height: 1.1 !important; }}"
    )


def make_checkbox_group(labels: list[str], active: list[int]) -> CheckboxGroup:
    return CheckboxGroup(labels=labels, active=active, styles=CHECKBOX_TEXT_STYLE)


def make_titled_panel(
    title: str,
    *children,
    width: int | None = None,
    styles: dict[str, str | None] | None = None,
):
    kwargs: dict[str, Any] = {"styles": styles} if styles is not None else {}
    if width is not None:
        kwargs["width"] = width
    return column(Div(text=f"<b>{title}</b>", styles=PANEL_TITLE_STYLE), *children, **kwargs)


def make_large_button(
    label: str,
    button_type: Literal["default", "primary", "success", "warning", "danger", "light"],
    width: int,
    height: int,
    button_text_css: str,
) -> Button:
    return Button(
        label=label,
        button_type=button_type,
        width=width,
        height=height,
        stylesheets=[InlineStyleSheet(css=button_text_css)],
    )


def js_add_visibility_bundle(js_parts: list[str], line_key: str, point_key: str, bar_key: str, condition: str):
    js_parts.append(f"{line_key}.visible = !!(({condition}) * (!use_bars));")
    js_parts.append(f"{point_key}.visible = !!(({condition}) * (!use_bars));")
    js_parts.append(f"{bar_key}.visible = !!(({condition}) * use_bars);")


def js_add_aggregate_from_sources(
    js_parts: list[str],
    out_src_key: str,
    source_terms: list[tuple[str | None, str]],
    var_tag: str,
    base_label: str,
    read_literal: str,
    reference_literal: str,
):
    pos_to_idx = f"pos_to_idx_{var_tag}"
    positions = f"positions_{var_tag}"
    y_count = f"y_count_{var_tag}"
    y_norm = f"y_norm_{var_tag}"
    sd_var = f"sd_{var_tag}"
    pv_var = f"pv_{var_tag}"
    ii_var = f"ii_{var_tag}"
    ji_var = f"ji_{var_tag}"
    order = f"order_{var_tag}"
    pos_sorted = f"pos_sorted_{var_tag}"
    y_count_sorted = f"y_count_sorted_{var_tag}"
    y_norm_sorted = f"y_norm_sorted_{var_tag}"

    js_parts.append(f"var {pos_to_idx} = {{}}; var {positions} = []; var {y_count} = []; var {y_norm} = [];")
    for cond_expr, src_key in source_terms:
        if cond_expr:
            js_parts.append(f"if ({cond_expr}) {{")
        js_parts.append(f"  var {sd_var} = {src_key}.data;")
        js_parts.append(f"  {sd_var}.position.forEach(function({pv_var}, {ii_var}) {{")
        js_parts.append(
            f"    if (!({pv_var} in {pos_to_idx})) {{ {pos_to_idx}[{pv_var}] = {positions}.length; {positions}.push({pv_var}); {y_count}.push(0); {y_norm}.push(0); }}"
        )
        js_parts.append(f"    var {ji_var} = {pos_to_idx}[{pv_var}];")
        js_parts.append(f"    {y_count}[{ji_var}] += {sd_var}.y_count[{ii_var}];")
        js_parts.append(f"    {y_norm}[{ji_var}] += {sd_var}.y_norm[{ii_var}];")
        js_parts.append("  });")
        if cond_expr:
            js_parts.append("}")

    js_parts.append(
        f"var {order} = {positions}.map(function(pv, i) {{ return {{pv: pv, i: i}}; }}).sort(function(a, b) {{ return a.pv - b.pv; }});"
    )
    js_parts.append(f"var {pos_sorted} = {order}.map(function(o) {{ return o.pv; }});")
    js_parts.append(f"var {y_count_sorted} = {order}.map(function(o) {{ return {y_count}[o.i]; }});")
    js_parts.append(f"var {y_norm_sorted} = {order}.map(function(o) {{ return {y_norm}[o.i]; }});")
    js_parts.append(f"{out_src_key}.data = {{")
    js_parts.append(f"  position: {pos_sorted},")
    js_parts.append(f"  y_count: {y_count_sorted},")
    js_parts.append(f"  y_norm: {y_norm_sorted},")
    js_parts.append(f"  y: (use_norm ? {y_norm_sorted} : {y_count_sorted}),")
    js_parts.append(f"  base_change: Array({pos_sorted}.length).fill('{base_label}'),")
    js_parts.append(f"  read_num: Array({pos_sorted}.length).fill({read_literal}),")
    js_parts.append(f"  reference_order: Array({pos_sorted}.length).fill({reference_literal}),")
    js_parts.append("};")
    js_parts.append(f"{out_src_key}.change.emit();")


def build_interaction_code(
    has_read_checkbox: bool,
    source_keys: list[str],
    aggregate_renderer_keys_individual: list[tuple[int, int, str, str, str, str]],
    aggregate_renderer_keys_combined: list[tuple[int, str, str, str, str]],
    aggregate_renderer_keys_all_reads_individual: list[tuple[int, str, str, str, str]],
    aggregate_renderer_keys_all_reads_combined: tuple[str, str, str, str],
    renderer_keys_individual: list[tuple[int, int, int, str, str, str]],
    renderer_keys_combined: list[tuple[int, int, str, str, str]],
    reads_available: list[int],
    reference_orders_available: list[Any],
) -> str:
    js_parts: list[str] = []
    js_parts.append("var use_norm = norm_cb.active.includes(0);")
    js_parts.append("var use_bars = plot_mode_cb.active.includes(0);")
    js_parts.append("var sum_mm = sum_mm_cb.active.includes(0);")
    if has_read_checkbox:
        js_parts.append("var sum_reads = sum_rd_cb.active.includes(0);")
    else:
        js_parts.append("var sum_reads = false;")
    js_parts.append("var mm_active = mismatch_cb.active;")
    js_parts.append("var ref_active = ref_order_cb.active;")
    if has_read_checkbox:
        js_parts.append("var rd_active = read_cb.active;")

    # Avoid && (becomes &amp;&amp;) and > (becomes &gt;); use * for boolean AND, !== for compare.
    js_parts.append("var use_combined = (ref_active.length !== 0) * (ref_active.length !== 1);")
    js_parts.append("if (plot.yaxis) { if (plot.yaxis.length !== 0) { plot.yaxis[0].axis_label = use_norm ? norm_label : raw_label; } }")

    # Switch all base sources between raw and normalized y arrays.
    for source_key in source_keys:
        js_parts.append(f"{source_key}.data.y = (use_norm ? {source_key}.data.y_norm : {source_key}.data.y_count).slice();")
        js_parts.append(f"{source_key}.change.emit();")

    # Build aggregate sources for individual-reference mode.
    for rn_idx, ref_idx, _, _, _, src_key_agg in aggregate_renderer_keys_individual:
        source_terms = []
        for bc_idx, bc in enumerate(MISMATCHES):
            src_key = f"src_r{reads_available[rn_idx]}_{bc.replace('>','_')}_rf{ref_idx}"
            source_terms.append((f"mm_active.includes({bc_idx})", src_key))
        js_add_aggregate_from_sources(
            js_parts=js_parts,
            out_src_key=src_key_agg,
            source_terms=source_terms,
            var_tag=f"agg_r{reads_available[rn_idx]}_rf{ref_idx}",
            base_label="SUM",
            read_literal=str(reads_available[rn_idx]),
            reference_literal=f"'{reference_orders_available[ref_idx]}'",
        )

    # Build aggregate sources for combined-reference mode.
    for rn_idx, _, _, _, src_key_agg_c in aggregate_renderer_keys_combined:
        source_terms = []
        for bc_idx, bc in enumerate(MISMATCHES):
            src_key_comb = f"src_r{reads_available[rn_idx]}_{bc.replace('>','_')}_comb"
            source_terms.append((f"mm_active.includes({bc_idx})", src_key_comb))
        js_add_aggregate_from_sources(
            js_parts=js_parts,
            out_src_key=src_key_agg_c,
            source_terms=source_terms,
            var_tag=f"agg_r{reads_available[rn_idx]}_comb",
            base_label="SUM",
            read_literal=str(reads_available[rn_idx]),
            reference_literal="'combined'",
        )

    # Build all-reads aggregate sources for individual-reference mode.
    for ref_idx, _, _, _, src_key_all in aggregate_renderer_keys_all_reads_individual:
        source_terms = []
        for rn_idx, rn in enumerate(reads_available):
            src_key_per_read = f"agg_src_r{rn}_rf{ref_idx}"
            if has_read_checkbox:
                source_terms.append((f"rd_active.includes({rn_idx})", src_key_per_read))
            else:
                source_terms.append((None, src_key_per_read))
        js_add_aggregate_from_sources(
            js_parts=js_parts,
            out_src_key=src_key_all,
            source_terms=source_terms,
            var_tag=f"agg_all_rf{ref_idx}",
            base_label="SUM_READS",
            read_literal="'ALL'",
            reference_literal=f"'{reference_orders_available[ref_idx]}'",
        )

    # Build all-reads aggregate source for combined-reference mode.
    key_l_all_c, key_c_all_c, key_b_all_c, src_key_all_c = aggregate_renderer_keys_all_reads_combined
    source_terms = []
    for rn_idx, rn in enumerate(reads_available):
        src_key_per_read_comb = f"agg_src_r{rn}_comb"
        if has_read_checkbox:
            source_terms.append((f"rd_active.includes({rn_idx})", src_key_per_read_comb))
        else:
            source_terms.append((None, src_key_per_read_comb))
    js_add_aggregate_from_sources(
        js_parts=js_parts,
        out_src_key=src_key_all_c,
        source_terms=source_terms,
        var_tag="agg_all_comb",
        base_label="SUM_READS",
        read_literal="'ALL'",
        reference_literal="'combined'",
    )

    # Individual mismatch renderers visibility.
    for rn_idx, bc_idx, ref_idx, key_l, key_c, key_b in renderer_keys_individual:
        if has_read_checkbox:
            cond = f"(!sum_mm) * (!use_combined) * (mm_active.includes({bc_idx})) * (rd_active.includes({rn_idx})) * (ref_active.includes({ref_idx}))"
        else:
            cond = f"(!sum_mm) * (!use_combined) * (mm_active.includes({bc_idx})) * (ref_active.includes({ref_idx}))"
        js_add_visibility_bundle(js_parts, key_l, key_c, key_b, cond)

    # Combined-reference mismatch renderers visibility.
    for rn_idx, bc_idx, key_l_comb, key_c_comb, key_b_comb in renderer_keys_combined:
        if has_read_checkbox:
            cond_comb = f"(!sum_mm) * (use_combined) * (mm_active.includes({bc_idx})) * (rd_active.includes({rn_idx}))"
        else:
            cond_comb = f"(!sum_mm) * (use_combined) * (mm_active.includes({bc_idx}))"
        js_add_visibility_bundle(js_parts, key_l_comb, key_c_comb, key_b_comb, cond_comb)

    # Aggregate renderers visibility in individual-reference mode.
    for rn_idx, ref_idx, key_l_agg, key_c_agg, key_b_agg, _ in aggregate_renderer_keys_individual:
        if has_read_checkbox:
            cond_agg = f"(sum_mm) * (!sum_reads) * (!use_combined) * (mm_active.length !== 0) * (rd_active.includes({rn_idx})) * (ref_active.includes({ref_idx}))"
        else:
            cond_agg = f"(sum_mm) * (!sum_reads) * (!use_combined) * (mm_active.length !== 0) * (ref_active.includes({ref_idx}))"
        js_add_visibility_bundle(js_parts, key_l_agg, key_c_agg, key_b_agg, cond_agg)

    # Aggregate renderers visibility in combined-reference mode.
    for rn_idx, key_l_agg_c, key_c_agg_c, key_b_agg_c, _ in aggregate_renderer_keys_combined:
        if has_read_checkbox:
            cond_agg_c = f"(sum_mm) * (!sum_reads) * (use_combined) * (mm_active.length !== 0) * (rd_active.includes({rn_idx}))"
        else:
            cond_agg_c = f"(sum_mm) * (!sum_reads) * (use_combined) * (mm_active.length !== 0)"
        js_add_visibility_bundle(js_parts, key_l_agg_c, key_c_agg_c, key_b_agg_c, cond_agg_c)

    # All-reads aggregate visibility in individual-reference mode.
    for ref_idx, key_l_all, key_c_all, key_b_all, _ in aggregate_renderer_keys_all_reads_individual:
        cond_all = f"(sum_mm) * (sum_reads) * (!use_combined) * (mm_active.length !== 0) * (ref_active.includes({ref_idx}))"
        js_add_visibility_bundle(js_parts, key_l_all, key_c_all, key_b_all, cond_all)

    # All-reads aggregate visibility in combined-reference mode.
    cond_all_c = "(sum_mm) * (sum_reads) * (use_combined) * (mm_active.length !== 0)"
    js_add_visibility_bundle(js_parts, key_l_all_c, key_c_all_c, key_b_all_c, cond_all_c)

    # Force DataRange1d to recompute bounds after y-values/visibility updates.
    js_parts.append("if (plot.y_range) {")
    js_parts.append("  if ('start' in plot.y_range) { plot.y_range.start = NaN; }")
    js_parts.append("  if ('end' in plot.y_range) { plot.y_range.end = NaN; }")
    js_parts.append("  plot.y_range.change.emit();")
    js_parts.append("}")

    return "\n".join(js_parts)


def prepare_plot_data(df: pd.DataFrame):
    is_insert = df["is_insert_mode"].iloc[0]
    y_raw_field = "count" if "count" in df.columns else "normalized_frequency"

    df_all = df.copy()
    df_all["source_base"] = df_all["base_change"].str.split(">").str[0]

    # Per read/reference_order/position/source_base denominator includes self-matches.
    denom_individual = df_all.groupby(
        ["read_num", "reference_order", "position", "source_base"], dropna=False
    )[y_raw_field].transform("sum")
    df_all["normalized_count"] = df_all[y_raw_field] / denom_individual.where(denom_individual != 0)
    df_all["normalized_count"] = df_all["normalized_count"].fillna(0.0)

    # Keep only mismatches (drop self-matches like A>A) for plotting.
    df_mm = df_all[df_all["base_change"].isin(MISMATCHES)].copy()
    reads_available = sorted(df_mm["read_num"].unique())
    reference_orders_available = sorted(df_mm["reference_order"].unique())

    return is_insert, y_raw_field, df_mm, df_all, reads_available, reference_orders_available


def build_base_sources(
    df_mm: pd.DataFrame,
    df_all: pd.DataFrame,
    reads_available,
    reference_orders_available,
    y_raw_field: str,
):
    sources = {}
    sources_combined = {}

    for rn in reads_available:
        sources[rn] = {}
        sources_combined[rn] = {}
        for bc in MISMATCHES:
            sources[rn][bc] = {}

            for ref_ord in reference_orders_available:
                sub = df_mm[
                    (df_mm["read_num"] == rn)
                    & (df_mm["base_change"] == bc)
                    & (df_mm["reference_order"] == ref_ord)
                ].sort_values("position")
                sources[rn][bc][ref_ord] = ColumnDataSource(data=dict(
                    position=sub["position"].tolist(),
                    y=sub[y_raw_field].tolist(),
                    y_count=sub[y_raw_field].tolist(),
                    y_norm=sub["normalized_count"].tolist(),
                    base_change=[bc] * len(sub),
                    read_num=[rn] * len(sub),
                    reference_order=[ref_ord] * len(sub),
                ))

            # Combined source: sum all reference_orders for this read and mismatch.
            sub_num = df_mm[(df_mm["read_num"] == rn) & (df_mm["base_change"] == bc)].groupby(
                "position", as_index=False
            )[y_raw_field].sum()

            sub_den = df_all[df_all["read_num"] == rn].copy()
            sub_den["source_base"] = sub_den["base_change"].str.split(">").str[0]
            sub_den = sub_den.groupby(["position", "source_base"], as_index=False)[y_raw_field].sum()
            source_base = bc.split(">")[0]
            sub_den = sub_den[sub_den["source_base"] == source_base][["position", y_raw_field]]
            sub_den = sub_den.rename(columns={y_raw_field: "denom"})

            sub_combined = sub_num.merge(sub_den, on="position", how="left")
            sub_combined["y_norm"] = sub_combined[y_raw_field] / sub_combined["denom"].where(
                sub_combined["denom"] != 0
            )
            sub_combined["y_norm"] = sub_combined["y_norm"].fillna(0.0)

            sources_combined[rn][bc] = ColumnDataSource(data=dict(
                position=sub_combined["position"].tolist(),
                y=sub_combined[y_raw_field].tolist(),
                y_count=sub_combined[y_raw_field].tolist(),
                y_norm=sub_combined["y_norm"].tolist(),
                base_change=[bc] * len(sub_combined),
                read_num=[rn] * len(sub_combined),
                reference_order=["combined"] * len(sub_combined),
            ))

    return sources, sources_combined


def build_render_collections(
    p,
    sources,
    sources_combined,
    reads_available,
    reference_orders_available,
):
    renderers: dict[Any, dict[str, dict[Any, tuple[Any, Any, Any]]]] = {}
    renderers_combined: dict[Any, dict[str, tuple[Any, Any, Any]]] = {}
    aggregate_sources: dict[Any, dict[Any, ColumnDataSource]] = {}
    aggregate_sources_combined: dict[Any, ColumnDataSource] = {}
    aggregate_sources_all_reads: dict[Any, ColumnDataSource] = {}
    aggregate_sources_all_reads_combined: ColumnDataSource | None = None
    aggregate_renderers: dict[Any, dict[Any, tuple[Any, Any, Any]]] = {}
    aggregate_renderers_combined: dict[Any, tuple[Any, Any, Any] | None] = {}
    aggregate_renderers_all_reads: dict[Any, tuple[Any, Any, Any]] = {}
    aggregate_renderers_all_reads_combined: tuple[Any, Any, Any] | None = None

    for rn in reads_available:
        renderers[rn] = {}
        renderers_combined[rn] = {}
        aggregate_sources[rn] = {}
        aggregate_sources_combined[rn] = ColumnDataSource(data=empty_plot_data())
        aggregate_renderers[rn] = {}
        aggregate_renderers_combined[rn] = None
        dash = READ_DASH_PATTERNS.get(rn, "dotted")
        marker = READ_MARKER_SHAPES.get(rn, "square")

        for bc in MISMATCHES:
            renderers[rn][bc] = {}

            # Individual renderers for each reference_order.
            for ref_ord in reference_orders_available:
                show_individual_on_load = len(reference_orders_available) == 1
                r, c, b = make_renderer_pair(
                    p=p,
                    source=sources[rn][bc][ref_ord],
                    color=PALETTE[bc],
                    dash=dash,
                    marker=marker,
                    line_alpha=0.8,
                    point_alpha=0.6,
                    visible=show_individual_on_load,
                )
                renderers[rn][bc][ref_ord] = (r, c, b)

            # Combined renderers (shown when multiple reference_orders selected).
            show_combined_on_load = len(reference_orders_available) > 1
            r_comb, c_comb, b_comb = make_renderer_pair(
                p=p,
                source=sources_combined[rn][bc],
                color=PALETTE[bc],
                dash=dash,
                marker=marker,
                line_alpha=0.8,
                point_alpha=0.6,
                visible=show_combined_on_load,
            )
            renderers_combined[rn][bc] = (r_comb, c_comb, b_comb)

        # Aggregate sources/renderers per selected reference_order.
        for ref_ord in reference_orders_available:
            aggregate_sources[rn][ref_ord] = ColumnDataSource(data=empty_plot_data())
            agg_line, agg_pts, agg_bar = make_renderer_pair(
                p=p,
                source=aggregate_sources[rn][ref_ord],
                color="#111111",
                dash=dash,
                marker=marker,
                line_alpha=0.9,
                point_alpha=0.7,
                visible=False,
            )
            aggregate_renderers[rn][ref_ord] = (agg_line, agg_pts, agg_bar)

        # Aggregate renderer for combined-reference mode.
        agg_line_comb, agg_pts_comb, agg_bar_comb = make_renderer_pair(
            p=p,
            source=aggregate_sources_combined[rn],
            color="#111111",
            dash=dash,
            marker=marker,
            line_alpha=0.9,
            point_alpha=0.7,
            visible=False,
        )
        aggregate_renderers_combined[rn] = (agg_line_comb, agg_pts_comb, agg_bar_comb)

    # Aggregate renderers for summing selected reads into one curve.
    for ref_ord in reference_orders_available:
        aggregate_sources_all_reads[ref_ord] = ColumnDataSource(data=empty_plot_data())
        agg_all_line, agg_all_pts, agg_all_bar = make_renderer_pair(
            p=p,
            source=aggregate_sources_all_reads[ref_ord],
            color="#8B0000",
            dash="solid",
            marker="diamond",
            line_alpha=0.95,
            point_alpha=0.75,
            visible=False,
        )
        aggregate_renderers_all_reads[ref_ord] = (agg_all_line, agg_all_pts, agg_all_bar)

    aggregate_sources_all_reads_combined = ColumnDataSource(data=empty_plot_data())
    agg_all_line_comb, agg_all_pts_comb, agg_all_bar_comb = make_renderer_pair(
        p=p,
        source=aggregate_sources_all_reads_combined,
        color="#8B0000",
        dash="solid",
        marker="diamond",
        line_alpha=0.95,
        point_alpha=0.75,
        visible=False,
    )
    aggregate_renderers_all_reads_combined = (agg_all_line_comb, agg_all_pts_comb, agg_all_bar_comb)

    return {
        "renderers": renderers,
        "renderers_combined": renderers_combined,
        "aggregate_sources": aggregate_sources,
        "aggregate_sources_combined": aggregate_sources_combined,
        "aggregate_sources_all_reads": aggregate_sources_all_reads,
        "aggregate_sources_all_reads_combined": aggregate_sources_all_reads_combined,
        "aggregate_renderers": aggregate_renderers,
        "aggregate_renderers_combined": aggregate_renderers_combined,
        "aggregate_renderers_all_reads": aggregate_renderers_all_reads,
        "aggregate_renderers_all_reads_combined": aggregate_renderers_all_reads_combined,
    }


def register_callback_artifacts(
    cb_args: dict[str, Any],
    reads_available,
    reference_orders_available,
    renderers,
    renderers_combined,
    aggregate_renderers,
    aggregate_renderers_combined,
    aggregate_renderers_all_reads,
    aggregate_renderers_all_reads_combined,
    sources,
    sources_combined,
    aggregate_sources,
    aggregate_sources_combined,
    aggregate_sources_all_reads,
    aggregate_sources_all_reads_combined,
):
    renderer_keys_individual = []  # list of (rn_idx, bc_idx, ref_idx, key_l, key_c, key_b)
    renderer_keys_combined = []    # list of (rn_idx, bc_idx, key_l, key_c, key_b)
    aggregate_renderer_keys_individual = []  # list of (rn_idx, ref_idx, key_l, key_c, key_b, src_key)
    aggregate_renderer_keys_combined = []    # list of (rn_idx, key_l, key_c, key_b, src_key)
    aggregate_renderer_keys_all_reads_individual = []  # list of (ref_idx, key_l, key_c, key_b, src_key)
    aggregate_renderer_keys_all_reads_combined = None  # tuple (key_l, key_c, key_b, src_key)
    source_keys = []

    # Register individual renderers.
    for rn_idx, rn in enumerate(reads_available):
        for bc_idx, bc in enumerate(MISMATCHES):
            for ref_idx, ref_ord in enumerate(reference_orders_available):
                line_r, circ_r, bar_r = renderers[rn][bc][ref_ord]
                key_l = f"r{rn}_{bc.replace('>','_')}_rf{ref_idx}_l"
                key_c = f"r{rn}_{bc.replace('>','_')}_rf{ref_idx}_c"
                key_b = f"r{rn}_{bc.replace('>','_')}_rf{ref_idx}_b"
                cb_args[key_l] = line_r
                cb_args[key_c] = circ_r
                cb_args[key_b] = bar_r
                renderer_keys_individual.append((rn_idx, bc_idx, ref_idx, key_l, key_c, key_b))
                source_key = f"src_r{rn}_{bc.replace('>','_')}_rf{ref_idx}"
                cb_args[source_key] = sources[rn][bc][ref_ord]
                source_keys.append(source_key)

    # Register combined renderers.
    for rn_idx, rn in enumerate(reads_available):
        for bc_idx, bc in enumerate(MISMATCHES):
            line_r_comb, circ_r_comb, bar_r_comb = renderers_combined[rn][bc]
            key_l_comb = f"r{rn}_{bc.replace('>','_')}_comb_l"
            key_c_comb = f"r{rn}_{bc.replace('>','_')}_comb_c"
            key_b_comb = f"r{rn}_{bc.replace('>','_')}_comb_b"
            cb_args[key_l_comb] = line_r_comb
            cb_args[key_c_comb] = circ_r_comb
            cb_args[key_b_comb] = bar_r_comb
            renderer_keys_combined.append((rn_idx, bc_idx, key_l_comb, key_c_comb, key_b_comb))
            source_key_comb = f"src_r{rn}_{bc.replace('>','_')}_comb"
            cb_args[source_key_comb] = sources_combined[rn][bc]
            source_keys.append(source_key_comb)

    # Register aggregate renderers/sources for individual reference_order mode.
    for rn_idx, rn in enumerate(reads_available):
        for ref_idx, ref_ord in enumerate(reference_orders_available):
            line_r_agg, circ_r_agg, bar_r_agg = aggregate_renderers[rn][ref_ord]
            key_l_agg = f"agg_r{rn}_rf{ref_idx}_l"
            key_c_agg = f"agg_r{rn}_rf{ref_idx}_c"
            key_b_agg = f"agg_r{rn}_rf{ref_idx}_b"
            src_key_agg = f"agg_src_r{rn}_rf{ref_idx}"
            cb_args[key_l_agg] = line_r_agg
            cb_args[key_c_agg] = circ_r_agg
            cb_args[key_b_agg] = bar_r_agg
            cb_args[src_key_agg] = aggregate_sources[rn][ref_ord]
            aggregate_renderer_keys_individual.append((rn_idx, ref_idx, key_l_agg, key_c_agg, key_b_agg, src_key_agg))

    # Register aggregate renderers/sources for combined-reference mode.
    for rn_idx, rn in enumerate(reads_available):
        line_r_agg_c, circ_r_agg_c, bar_r_agg_c = aggregate_renderers_combined[rn]
        key_l_agg_c = f"agg_r{rn}_comb_l"
        key_c_agg_c = f"agg_r{rn}_comb_c"
        key_b_agg_c = f"agg_r{rn}_comb_b"
        src_key_agg_c = f"agg_src_r{rn}_comb"
        cb_args[key_l_agg_c] = line_r_agg_c
        cb_args[key_c_agg_c] = circ_r_agg_c
        cb_args[key_b_agg_c] = bar_r_agg_c
        cb_args[src_key_agg_c] = aggregate_sources_combined[rn]
        aggregate_renderer_keys_combined.append((rn_idx, key_l_agg_c, key_c_agg_c, key_b_agg_c, src_key_agg_c))

    # Register all-reads aggregate renderers/sources for individual reference_order mode.
    for ref_idx, ref_ord in enumerate(reference_orders_available):
        line_r_all, circ_r_all, bar_r_all = aggregate_renderers_all_reads[ref_ord]
        key_l_all = f"agg_all_rf{ref_idx}_l"
        key_c_all = f"agg_all_rf{ref_idx}_c"
        key_b_all = f"agg_all_rf{ref_idx}_b"
        src_key_all = f"agg_all_src_rf{ref_idx}"
        cb_args[key_l_all] = line_r_all
        cb_args[key_c_all] = circ_r_all
        cb_args[key_b_all] = bar_r_all
        cb_args[src_key_all] = aggregate_sources_all_reads[ref_ord]
        aggregate_renderer_keys_all_reads_individual.append((ref_idx, key_l_all, key_c_all, key_b_all, src_key_all))

    # Register all-reads aggregate renderer/source for combined-reference mode.
    line_r_all_c, circ_r_all_c, bar_r_all_c = aggregate_renderers_all_reads_combined
    key_l_all_c = "agg_all_comb_l"
    key_c_all_c = "agg_all_comb_c"
    key_b_all_c = "agg_all_comb_b"
    src_key_all_c = "agg_all_src_comb"
    cb_args[key_l_all_c] = line_r_all_c
    cb_args[key_c_all_c] = circ_r_all_c
    cb_args[key_b_all_c] = bar_r_all_c
    cb_args[src_key_all_c] = aggregate_sources_all_reads_combined
    aggregate_renderer_keys_all_reads_combined = (key_l_all_c, key_c_all_c, key_b_all_c, src_key_all_c)

    return {
        "renderer_keys_individual": renderer_keys_individual,
        "renderer_keys_combined": renderer_keys_combined,
        "aggregate_renderer_keys_individual": aggregate_renderer_keys_individual,
        "aggregate_renderer_keys_combined": aggregate_renderer_keys_combined,
        "aggregate_renderer_keys_all_reads_individual": aggregate_renderer_keys_all_reads_individual,
        "aggregate_renderer_keys_all_reads_combined": aggregate_renderer_keys_all_reads_combined,
        "source_keys": source_keys,
    }


def build_controls_layout(
    has_read_checkbox: bool,
    size_panel,
    plot_mode_panel,
    sum_panel,
    normalization_panel,
    reference_order_panel,
    mismatch_panel,
    controls_style: dict[str, str],
    panel_text_style: dict[str, str | None],
    read_checkbox=None,
    sum_reads_panel=None,
    interaction_cb=None,
):
    agg_children = [sum_panel]
    if has_read_checkbox and sum_reads_panel is not None:
        agg_children.append(sum_reads_panel)
    agg_children.append(normalization_panel)

    agg_stack = column(
        *agg_children,
        spacing=12,
        width=340,
        css_classes=["tm-agg-stack"],
    )

    if has_read_checkbox and read_checkbox is not None:
        read_panel = make_titled_panel(
            "Reads",
            read_checkbox,
            width=200,
            styles=panel_text_style,
        )
        if interaction_cb is not None:
            read_checkbox.js_on_change("active", interaction_cb)
        return row(
            size_panel,
            Spacer(width=20),
            plot_mode_panel,
            Spacer(width=20),
            agg_stack,
            Spacer(width=30),
            read_panel,
            Spacer(width=20),
            reference_order_panel,
            Spacer(width=20),
            mismatch_panel,
            sizing_mode="stretch_width",
            styles=controls_style,
        )

    return row(
        size_panel,
        Spacer(width=20),
        plot_mode_panel,
        Spacer(width=20),
        agg_stack,
        Spacer(width=20),
        reference_order_panel,
        Spacer(width=20),
        mismatch_panel,
        sizing_mode="stretch_width",
        styles=controls_style,
    )


def create_plot_figure(is_insert: bool, y_raw_field: str):
    x_label = "Fragment position" if is_insert else "Read position"
    p: Any = figure(
        title="Tasmanian-mismatch profile",
        x_axis_label=x_label,
        y_axis_label=y_raw_field.replace("_", " ").title(),
        sizing_mode="stretch_both",
        min_width=760,
        min_height=420,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    p.x_range.only_visible = True
    p.y_range.only_visible = True

    p.title.text_font_size = "18pt"
    p.xaxis.axis_label_text_font_size = "16pt"
    p.yaxis.axis_label_text_font_size = "16pt"
    p.xaxis.major_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = "14pt"

    hover = HoverTool(tooltips=[
        ("Position", "@position"),
        ("Displayed Value", "@y{0,0.000000}"),
        ("Count", "@y_count{0,0.000000}"),
        ("Normalized", "@y_norm{0,0.000000}"),
        ("Mismatch", "@base_change"),
        ("Read", "@read_num"),
    ])
    p.add_tools(hover)
    return p


def write_plot_html(layout, output_path: str | None = None, extra_js: str = ""):
    html = file_html(
        layout,
        resources=Resources(mode="cdn"),
        title="Tasmanian-mismatch",
    )
    style_inject = (
        "<style>"
        ":root {"
        "  --tm-bg-1: #f6f8fb;"
        "  --tm-bg-2: #edf2f8;"
        "  --tm-card: rgba(255, 255, 255, 0.86);"
        "  --tm-card-border: rgba(19, 31, 50, 0.10);"
        "  --tm-title: #1a2a3f;"
        "  --tm-text: #2f3f58;"
        "  --tm-accent: #0f6bd6;"
        "}"
        "html, body {"
        "  margin:0; padding:0; width:100%; min-height:100%;"
        "  background: linear-gradient(120deg, var(--tm-bg-1), var(--tm-bg-2));"
        "  color: var(--tm-text);"
        "  font-family: 'Segoe UI', 'IBM Plex Sans', 'Helvetica Neue', sans-serif;"
        "}"
        ".bk-root {"
        "  padding: 10px 14px 18px 14px;"
        "  border-radius: 14px;"
        "  background: var(--tm-card);"
        "  box-shadow: 0 12px 30px rgba(12, 24, 40, 0.08);"
        "  border: 1px solid var(--tm-card-border);"
        "}"
        ".bk-root .bk-btn {"
        "  border-radius: 20px !important;"
        "  border: 1px solid rgba(28, 43, 66, 0.14);"
        "  background: #ffffff;"
        "  color: #21324a;"
        "  font-weight: 700 !important;"
        "  font-size: 28px !important;"
        "  line-height: 1.2 !important;"
        "  letter-spacing: 0.1px;"
        "  transition: transform 140ms ease, box-shadow 140ms ease;"
        "  padding: 2px 10px !important;"
        "}"
        ".bk-root .bk-btn span, .bk-root .bk-btn-text {"
        "  font-size: 28px !important;"
        "  font-weight: 700 !important;"
        "}"
        ".bk-root .bk-btn:hover {"
        "  transform: translateY(-1px);"
        "  box-shadow: 0 6px 16px rgba(12, 24, 40, 0.14);"
        "}"
        ".bk-root .bk-input-group label {"
        "  color: var(--tm-text);"
        "  font-weight: 500;"
        "  letter-spacing: 0.12px;"
        "}"
        ".bk-root .tm-agg-stack .bk-input-group label {"
        "  display: block;"
        "  max-width: 300px;"
        "  white-space: normal;"
        "  overflow-wrap: anywhere;"
        "  line-height: 1.35;"
        "}"
        ".bk-root input[type='checkbox'] { accent-color: var(--tm-accent); }"
        ".bk-root .bk-slider-title {"
        "  color: var(--tm-title);"
        "  font-weight: 650;"
        "  letter-spacing: 0.2px;"
        "}"
        ".bk-root .tm-resizable-plot {"
        "  position: relative;"
        "  width: 1280px;"
        "  height: 620px;"
        "  min-width: 760px;"
        "  min-height: 420px;"
        "  max-width: 100%;"
        "  overflow: auto;"
        "  border: 1px solid rgba(19,31,50,0.16);"
        "  border-radius: 12px;"
        "  background: rgba(255,255,255,0.82);"
        "  box-shadow: inset 0 0 0 1px rgba(255,255,255,0.3);"
        "  box-sizing: border-box;"
        "}"
        ".bk-root .tm-resizable-plot::after {"
        "  content: '↙';"
        "  position: absolute;"
        "  right: 6px;"
        "  bottom: 6px;"
        "  width: 30px;"
        "  height: 30px;"
        "  display: flex;"
        "  align-items: center;"
        "  justify-content: center;"
        "  cursor: se-resize;"
        "  z-index: 50;"
        "  opacity: 0.6;"
        "  transition: opacity 120ms, transform 120ms;"
        "  font-size: 22px;"
        "  color: #30445f;"
        "  font-weight: bold;"
        "  user-select: none;"
        "  pointer-events: auto;"
        "  background: rgba(255,255,255,0.85);"
        "  border: 1px solid rgba(19,31,50,0.16);"
        "  border-radius: 6px;"
        "  box-shadow: 0 2px 6px rgba(0,0,0,0.12);"
        "}"
        ".bk-root .tm-resizable-plot::after:hover {"
        "  opacity: 0.95;"
        "  transform: scale(1.08);"
        "}"
        "@media (max-width: 960px) {"
        "  .bk-root { padding: 8px; border-radius: 10px; }"
        "  .bk-root .tm-resizable-plot { width: 100%; height: 520px; }"
        "}"
        "</style>"
    )
    html = html.replace("</head>", style_inject + "</head>")
    html = html.replace("</body>", '<div style="height:80px"></div></body>')
    if extra_js:
        html = html.replace("</body>", extra_js + "</body>")

    if output_path:
        Path(output_path).write_text(html)
        print(f"Saved to {output_path}")
        return

    import tempfile
    import webbrowser

    tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False, mode="w")
    tmp.write(html)
    tmp.close()
    print(f"Opening {tmp.name}")
    webbrowser.open(f"file://{tmp.name}")


def build_plot(df: pd.DataFrame, output_path: str | None = None):
    is_insert, y_raw_field, df, df_all, reads_available, reference_orders_available = prepare_plot_data(df)
    sources, sources_combined = build_base_sources(
        df_mm=df,
        df_all=df_all,
        reads_available=reads_available,
        reference_orders_available=reference_orders_available,
        y_raw_field=y_raw_field,
    )
    p = create_plot_figure(is_insert=is_insert, y_raw_field=y_raw_field)

    # ---- Render lines ------------------------------------------------------
    # Read 1 = solid, Read 2 = dashed (so overlapping lines are distinguishable)
    render_data = build_render_collections(
        p=p,
        sources=sources,
        sources_combined=sources_combined,
        reads_available=reads_available,
        reference_orders_available=reference_orders_available,
    )
    renderers = render_data["renderers"]
    renderers_combined = render_data["renderers_combined"]
    aggregate_sources = render_data["aggregate_sources"]
    aggregate_sources_combined = render_data["aggregate_sources_combined"]
    aggregate_sources_all_reads = render_data["aggregate_sources_all_reads"]
    aggregate_sources_all_reads_combined = render_data["aggregate_sources_all_reads_combined"]
    aggregate_renderers = render_data["aggregate_renderers"]
    aggregate_renderers_combined = render_data["aggregate_renderers_combined"]
    aggregate_renderers_all_reads = render_data["aggregate_renderers_all_reads"]
    aggregate_renderers_all_reads_combined = render_data["aggregate_renderers_all_reads_combined"]

    # ---- Mismatch checkbox group with Select All / Deselect All -----------
    mismatch_checkbox = make_checkbox_group(
        labels=MISMATCHES,
        active=list(range(len(MISMATCHES))),
    )
    button_text_css = build_button_text_css(font_px=28)
    select_all_btn = make_large_button("Select All", "primary", 180, 40, button_text_css)
    deselect_all_btn = make_large_button("Deselect All", "default", 190, 40, button_text_css)
    # JS to select / deselect all mismatches (no arrow functions - Bokeh escapes > to &gt;)
    select_all_btn.js_on_click(CustomJS(args={"cb": mismatch_checkbox}, code=
        "cb.active = Array.from({length: cb.labels.length}, function(_, i) { return i; });"))
    deselect_all_btn.js_on_click(CustomJS(args={"cb": mismatch_checkbox}, code=
        "cb.active = [];"))

    # ---- Build unified interaction callback --------------------------------
    # Reference-order selection can switch between individual and combined-reference mode.
    # Optional mismatch-sum mode collapses selected mismatch types into one curve.

    # Create reference_order checkbox
    reference_order_checkbox = make_checkbox_group(
        labels=[str(ro) for ro in reference_orders_available],
        active=list(range(len(reference_orders_available))),
    )

    normalize_checkbox = make_checkbox_group(labels=["Normalize counts"], active=[])
    plot_mode_checkbox = make_checkbox_group(labels=["Bar mode"], active=[])

    sum_mismatch_checkbox = make_checkbox_group(
        labels=["Sum shown mismatches into one curve"],
        active=[],
    )

    # Collect all renderer references into a flat JS-accessible structure.
    cb_args = {
        "mismatch_cb": mismatch_checkbox,
        "ref_order_cb": reference_order_checkbox,
        "norm_cb": normalize_checkbox,
        "plot_mode_cb": plot_mode_checkbox,
        "sum_mm_cb": sum_mismatch_checkbox,
        "plot": p,
        "raw_label": y_raw_field.replace("_", " ").title(),
        "norm_label": "Normalized Count",
    }

    has_read_checkbox = not is_insert and len(reads_available) >= 2

    if has_read_checkbox:
        read_checkbox = make_checkbox_group(
            labels=[f"Read {rn}" for rn in reads_available],
            active=list(range(len(reads_available))),
        )
        sum_reads_checkbox = make_checkbox_group(
            labels=["Sum selected reads into one curve"],
            active=[],
        )
        cb_args["read_cb"] = read_checkbox
        cb_args["sum_rd_cb"] = sum_reads_checkbox

    callback_artifacts = register_callback_artifacts(
        cb_args=cb_args,
        reads_available=reads_available,
        reference_orders_available=reference_orders_available,
        renderers=renderers,
        renderers_combined=renderers_combined,
        aggregate_renderers=aggregate_renderers,
        aggregate_renderers_combined=aggregate_renderers_combined,
        aggregate_renderers_all_reads=aggregate_renderers_all_reads,
        aggregate_renderers_all_reads_combined=aggregate_renderers_all_reads_combined,
        sources=sources,
        sources_combined=sources_combined,
        aggregate_sources=aggregate_sources,
        aggregate_sources_combined=aggregate_sources_combined,
        aggregate_sources_all_reads=aggregate_sources_all_reads,
        aggregate_sources_all_reads_combined=aggregate_sources_all_reads_combined,
    )
    renderer_keys_individual = callback_artifacts["renderer_keys_individual"]
    renderer_keys_combined = callback_artifacts["renderer_keys_combined"]
    aggregate_renderer_keys_individual = callback_artifacts["aggregate_renderer_keys_individual"]
    aggregate_renderer_keys_combined = callback_artifacts["aggregate_renderer_keys_combined"]
    aggregate_renderer_keys_all_reads_individual = callback_artifacts["aggregate_renderer_keys_all_reads_individual"]
    aggregate_renderer_keys_all_reads_combined = callback_artifacts["aggregate_renderer_keys_all_reads_combined"]
    source_keys = callback_artifacts["source_keys"]

    # Build one callback - NO < > or => anywhere: Bokeh HTML-escapes them to &lt; &gt; &gt;
    # breaking the JS. Use !== for comparisons, forEach for loops, function() for lambdas.
    interaction_code = build_interaction_code(
        has_read_checkbox=has_read_checkbox,
        source_keys=source_keys,
        aggregate_renderer_keys_individual=aggregate_renderer_keys_individual,
        aggregate_renderer_keys_combined=aggregate_renderer_keys_combined,
        aggregate_renderer_keys_all_reads_individual=aggregate_renderer_keys_all_reads_individual,
        aggregate_renderer_keys_all_reads_combined=aggregate_renderer_keys_all_reads_combined,
        renderer_keys_individual=renderer_keys_individual,
        renderer_keys_combined=renderer_keys_combined,
        reads_available=reads_available,
        reference_orders_available=reference_orders_available,
    )
    interaction_cb = CustomJS(args=cb_args, code=interaction_code)

    # Wire all checkboxes to the same callback.
    mismatch_checkbox.js_on_change("active", interaction_cb)
    reference_order_checkbox.js_on_change("active", interaction_cb)
    normalize_checkbox.js_on_change("active", interaction_cb)
    plot_mode_checkbox.js_on_change("active", interaction_cb)
    sum_mismatch_checkbox.js_on_change("active", interaction_cb)
    if has_read_checkbox:
        sum_reads_checkbox.js_on_change("active", interaction_cb)
    # Select/deselect all — no arrow functions (=> uses >)
    select_all_btn.js_on_click(CustomJS(args=cb_args, code=
        "mismatch_cb.active = Array.from({length: mismatch_cb.labels.length}, function(_, i) { return i; });\n"
        + interaction_code))
    deselect_all_btn.js_on_click(CustomJS(args=cb_args, code=
        "mismatch_cb.active = [];\n" + interaction_code))

    # ---- Color legend (non-interactive, just for reference) ----------------
    legend_items = []
    for bc in MISMATCHES:
        rn = reads_available[0]
        line_r, _, _ = renderers_combined[rn][bc]
        legend_items.append(LegendItem(label=bc, renderers=[line_r]))
    legend = Legend(items=legend_items, click_policy="none", location="top_right")
    legend.label_text_font_size = "13pt"
    p.add_layout(legend, "right")

    # ---- Layout -----------------------------------------------------------
    panel_text_style: dict[str, str | None] = {"font-size": "14pt"}

    mismatch_panel = make_titled_panel(
        "Mismatch types",
        row(select_all_btn, deselect_all_btn),
        mismatch_checkbox,
        width=320,
        styles=panel_text_style,
    )

    reference_order_panel = make_titled_panel(
        "Reference Order",
        reference_order_checkbox,
        width=280,
        styles=panel_text_style,
    )

    normalization_panel = make_titled_panel(
        "Y Scale",
        normalize_checkbox,
        width=220,
        styles=panel_text_style,
    )

    plot_mode_panel = make_titled_panel(
        "Plot Style",
        plot_mode_checkbox,
        width=220,
        styles=panel_text_style,
    )

    sum_panel = make_titled_panel(
        "Mismatch Agg.",
        sum_mismatch_checkbox,
        width=200,
        styles=panel_text_style,
    )

    plot_container = column(
        p,
        css_classes=["tm-resizable-plot"],
        sizing_mode="fixed",
        width=1280,
        height=620,
    )

    size_panel = make_titled_panel(
        "Plot Size",
        Div(
            text=(
                "<div id='tm-plot-grip' style='width:54px;height:54px;display:flex;align-items:center;"
                "justify-content:center;cursor:nwse-resize;font-size:34px;color:#30445f;font-weight:800;"
                "background:rgba(255,255,255,0.95);border:2px solid rgba(19,31,50,0.25);"
                "border-radius:12px;box-shadow:0 2px 10px rgba(0,0,0,0.10);user-select:none;"
                "pointer-events:auto;touch-action:none;'>"
                "&#x2198;"
                "</div>"
                "<div style='font-size:15px;line-height:1.35;color:#354861;margin-top:6px;'>"
                "Drag to resize plot"
                "</div>"
            ),
            css_classes=["tm-resize-grip-panel"],
        ),
        width=220,
    )

    if has_read_checkbox:
        sum_reads_panel = make_titled_panel(
            "Read Agg.",
            sum_reads_checkbox,
            width=200,
            styles=panel_text_style,
        )

    # Horizontal control panel - spread widgets across width
    CONTROLS_STYLE = {
        "gap": "32px",
        "padding": "16px 20px",
        "border-bottom": "none",
        "background": "rgba(240,242,245,0.4)",
        "overflow-x": "auto",
    }

    controls = build_controls_layout(
        has_read_checkbox=has_read_checkbox,
        size_panel=size_panel,
        plot_mode_panel=plot_mode_panel,
        sum_panel=sum_panel,
        normalization_panel=normalization_panel,
        reference_order_panel=reference_order_panel,
        mismatch_panel=mismatch_panel,
        controls_style=CONTROLS_STYLE,
        panel_text_style=panel_text_style,
        read_checkbox=read_checkbox if has_read_checkbox else None,
        sum_reads_panel=sum_reads_panel if has_read_checkbox else None,
        interaction_cb=interaction_cb,
    )

    # Keep controls above the plot to guarantee visibility in standalone HTML.
    layout = column(controls, plot_container, sizing_mode="stretch_width")

    # Inject a deterministic in-plot resize grip and keep Bokeh model dimensions
    # synchronized with container size changes.
    container_id = plot_container.id
    resize_js = f"""
<script>
(function() {{
    var CID = '{container_id}';
    var retryCount = 0;
    var maxRetries = 40;

    function deepFind(root, selector) {{
        var found = root.querySelector(selector);
        if (found) return found;
        var all = root.querySelectorAll('*');
        for (var i = 0; i < all.length; i++) {{
            if (all[i].shadowRoot) {{
                var r = deepFind(all[i].shadowRoot, selector);
                if (r) return r;
            }}
        }}
        return null;
    }}

    function findGrip() {{
        return deepFind(document, '#tm-plot-grip');
    }}

    // Bokeh buttons render in shadow roots. Force larger label text for specific controls.
    function enforceButtonFonts() {{
        var targetLabels = {{
            'Select All': true,
            'Deselect All': true,
        }};
        var changed = 0;
        var all = document.querySelectorAll('*');
        for (var i = 0; i < all.length; i++) {{
            var host = all[i];
            if (!host.shadowRoot) continue;
            var btns = host.shadowRoot.querySelectorAll('button, .bk-btn');
            for (var j = 0; j < btns.length; j++) {{
                var b = btns[j];
                var txt = (b.textContent || '').trim();
                if (!targetLabels[txt]) continue;
                b.style.fontSize = '28px';
                b.style.fontWeight = '800';
                b.style.lineHeight = '1.05';
                changed++;
            }}
        }}
        if (changed > 0) console.log('[TM] Button font overrides applied to', changed, 'buttons');
    }}
    // Run multiple times because Bokeh mounts widgets asynchronously.
    for (var pass = 0; pass < 12; pass++) {{
        (function(delay) {{ setTimeout(enforceButtonFonts, delay); }})(pass * 200);
    }}

    function getModel() {{
        try {{
            if (typeof Bokeh === 'undefined') return null;
            var docs = Bokeh.documents;
            if (!docs || !docs.length) return null;
            for (var i = 0; i < docs.length; i++) {{
                var m = docs[i].get_model_by_id(CID);
                if (m) return m;
            }}
        }} catch(e) {{}}
        return null;
    }}

    function wireGrip() {{
        retryCount++;
        var model = getModel();
        if (!model) {{
            if (retryCount < maxRetries) {{
                setTimeout(wireGrip, 200);
                return;
            }} else {{
                console.error('[TM] Could not find Bokeh model', CID, 'after', retryCount, 'retries');
                return;
            }}
        }}
        
        var grip = findGrip();
        if (!grip) {{
            if (retryCount < maxRetries) {{
                setTimeout(wireGrip, 200);
                return;
            }} else {{
                console.error('[TM] Could not find grip after', retryCount, 'retries.');
                return;
            }}
        }}
        
        console.log('[TM] Grip found, wiring drag. Initial plot size:', model.width, 'x', model.height);

        grip.style.opacity = '0.8';
        grip.onmouseenter = function() {{ grip.style.opacity = '1'; grip.style.transform = 'scale(1.15)'; }};
        grip.onmouseleave = function() {{ grip.style.opacity = '0.8'; grip.style.transform = ''; }};

        grip.addEventListener('mousedown', function(evt) {{
            console.log('[TM] Mousedown on grip');
            var startX = evt.clientX;
            var startY = evt.clientY;
            var startW = model.width  || 1280;
            var startH = model.height || 620;
            console.log('[TM] Drag start:', startW, 'x', startH);
            evt.preventDefault();
            evt.stopPropagation();

            function onMove(e) {{
                var newW = Math.max(600, startW + (e.clientX - startX));
                var newH = Math.max(350, startH + (e.clientY - startY));
                model.width  = newW;
                model.height = newH;
                model.change.emit();
            }}
            function onUp() {{
                console.log('[TM] Drag end:', model.width, 'x', model.height);
                document.removeEventListener('mousemove', onMove);
                document.removeEventListener('mouseup', onUp);
            }}
            document.addEventListener('mousemove', onMove);
            document.addEventListener('mouseup', onUp);
        }});
        console.log('[TM] Drag listener attached to grip element');
    }}
    wireGrip();
}})();
</script>"""

    write_plot_html(layout=layout, output_path=output_path, extra_js=resize_js)


def main():
    parser = argparse.ArgumentParser(description="Visualize tasmanian-mismatch output with Bokeh")
    parser.add_argument("input", help="TSV output from tasmanian-mismatch")
    parser.add_argument("-o", "--output", help="Save standalone HTML instead of opening browser")
    args = parser.parse_args()

    df = load_data(args.input)
    build_plot(df, args.output)


if __name__ == "__main__":
    main()
