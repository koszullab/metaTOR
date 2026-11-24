#!/usr/bin/env python3
# coding: utf-8

"""The report of a run: what the map of its clusters holds, digested, for a
person to read and for MultiQC to aggregate.

metator report writes it in an output directory from the map of the clusters
of a run alone (see collect), and so do metator export and pipeline, in its
step export:

- report.txt (REPORT_TXT): a text file, in sections, a key and its value per
  line, then tables, a row per line; keys, columns and units are fixed by
  FORMAT_VERSION, so that a script or a MultiQC module parses it (see
  write_txt and read_txt);
- report_mqc.yaml (REPORT_MQC): MultiQC's custom content, a row of its
  General Statistics table, which ``multiqc`` finds and aggregates over the
  runs (see write_mqc);
- report.html (REPORT_HTML): a page to read the run with, from its Jinja2
  template (templates/report.html.j2), its plots drawn by plotly.js, which the
  page loads from a CDN (PLOTLY_JS) when it is opened (see write_html).

The quality classes of the final clusters (HQ, MQ, LQ, Contaminated, Other:
see metator.validation.bin_quality) are those of the values of miComplete,
which the map holds at full precision, as the log of the run counts them.
"""

import itertools
import math
import os
import re
import time
from dataclasses import dataclass, field
from os.path import join

import numpy as np
import pandas as pd

import metator.clusters as mtcl
import metator.figures as mtf
import metator.host as mth
import metator.metacool as mtmc
import metator.outputs as mto
import metator.validation as mtv
from metator.version import __version__

#: The version of the format of report.txt and report_mqc.yaml: a key, a
#: column or a unit that changes, or goes, changes it; one that comes does
#: not.
FORMAT_VERSION = 1

REPORT_TXT = "report.txt"
REPORT_HTML = "report.html"
REPORT_MQC = "report_mqc.yaml"
REPORT_FILES = (REPORT_TXT, REPORT_HTML, REPORT_MQC)

#: plotly.js, which report.html loads from its CDN when it is opened: a
#: pinned version, which draws the plots the same over time.
PLOTLY_JS = "https://cdn.plot.ly/plotly-2.35.2.min.js"

#: The quality classes of the final clusters, in their order (see
#: metator.validation.bin_quality), and their colours in report.html.
QUALITIES = ("HQ", "MQ", "LQ", "Contaminated", "Other")
QUALITY_COLOURS = {"HQ": "#1a9850", "MQ": "#91cf60", "LQ": "#fdae61", "Contaminated": "#d73027", "Other": "#bdbdbd"}

#: The steps of a run that the map records, in their order, under the metator
#: key of its metadata.
STEPS = ("network", "partition", "validation", "mge", "host")

#: The most contigs that the plot of the GC content and the coverage of the
#: contigs of report.html shows: those of the MAGs and of the MGE clusters
#: first, the longest first, then the others.
MAX_PLOTTED_CONTIGS = 20000

#: The most MAGs that the heatmap of their contacts of report.html shows: the
#: largest, by quality class.
MAX_HEATMAP_MAGS = 150

#: The columns of the table of the MAGs of report.txt.
MAG_COLUMNS = (
    "name",
    "quality",
    "length_bp",
    "contigs",
    "completeness_pct",
    "contamination_pct",
    "gc_pct",
    "n50_bp",
    "hic_abundance",
    "sg_abundance",
    "intra_genome_ratio",
    "n_mges",
    "mges",
)

#: The columns of the General Statistics table of MultiQC (report_mqc.yaml):
#: the key of report.txt, its section, and how MultiQC shows it.
MQC_COLUMNS = (
    (
        "contigs_bp",
        "contigs",
        {"title": "Assembly", "description": "Length of the contigs of the map", "suffix": " bp", "format": "{:,.0f}"},
    ),
    (
        "contig_n50_bp",
        "contigs",
        {"title": "Contig N50", "description": "N50 of the contigs of the map", "suffix": " bp", "format": "{:,.0f}"},
    ),
    (
        "trans_contig_ratio",
        "contacts",
        {
            "title": "3D ratio",
            "description": "Share of the contacts between two contigs",
            "min": 0,
            "max": 1,
            "format": "{:,.3f}",
        },
    ),
    (
        "in_mags_pct",
        "contigs",
        {
            "title": "Binned",
            "description": "Share of the length of the contigs in a MAG",
            "suffix": "%",
            "min": 0,
            "max": 100,
            "format": "{:,.1f}",
        },
    ),
    (
        "hq_mags",
        "clusters",
        {
            "title": "HQ MAGs",
            "description": "Final clusters of completeness >= 90% and contamination <= 5%",
            "format": "{:,.0f}",
            "scale": "Greens",
        },
    ),
    (
        "mq_mags",
        "clusters",
        {
            "title": "MQ MAGs",
            "description": "Final clusters of completeness >= 70% and contamination <= 10%",
            "format": "{:,.0f}",
            "scale": "Greens",
        },
    ),
    (
        "lq_mags",
        "clusters",
        {
            "title": "LQ MAGs",
            "description": "Final clusters of completeness >= 50% and contamination <= 10%",
            "format": "{:,.0f}",
        },
    ),
    (
        "contaminated_mags",
        "clusters",
        {
            "title": "Contaminated",
            "description": "Final clusters of completeness >= 50% and contamination > 10%",
            "format": "{:,.0f}",
            "scale": "Reds",
        },
    ),
    (
        "mge_clusters",
        "mges",
        {"title": "MGE clusters", "description": "Clusters of the MGE contigs (metator mge)", "format": "{:,.0f}"},
    ),
    (
        "mge_clusters_with_host",
        "hosts",
        {"title": "Hosted MGEs", "description": "MGE clusters with at least one host MAG (metator host)", "format": "{:,.0f}"},
    ),
)


@dataclass
class Report:
    """The report of a run (see collect).

    Attributes
    ----------
    sections : dict
        The sections of report.txt, by name, in order: each its keys and
        their values, in order (int, float, str; None when missing).
    tables : dict
        The tables of report.txt, by name, in order (pandas.DataFrame).
    plots : dict
        What the plots of report.html are drawn from, by plot.
    """

    sections: dict = field(default_factory=dict)
    tables: dict = field(default_factory=dict)
    plots: dict = field(default_factory=dict)


def sample_name(metacool):
    """The name of the run in its report, which MultiQC names its row with: that
    of the file of its map without .metacool, or of its entry for an entry of a
    collection (file.mmetacool::entries/NAME)."""
    if mtmc.is_entry(metacool):
        name = str(metacool).split("::", 1)[1].strip().strip("/").rsplit("/", 1)[-1]
    else:
        name = os.path.basename(mtmc.map_file(metacool))
    return re.sub(r"\.m?metacool$", "", name) or name


def _n50(lengths):
    """The N50 and the L50 of lengths (0 and 0 for none)."""
    ordered = np.sort(np.asarray(lengths, dtype=np.int64))[::-1]
    if not len(ordered):
        return 0, 0
    reached = np.cumsum(ordered) >= ordered.sum() / 2
    position = int(np.argmax(reached))
    return int(ordered[position]), position + 1


def _number(value):
    """A value of the map as a Python number: None for a missing one (NaN,
    None, an empty text)."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _names(value):
    """A list of names stored as a text, separated by commas (the hosts of an
    MGE cluster, the MGE clusters of a MAG): none for a missing one."""
    if not isinstance(value, str):
        return []
    return [name for name in value.split(mth.SEPARATOR) if name and name != "nan"]


def _percent(part, whole):
    """part of whole, in percent: None for an empty whole."""
    return None if not whole else 100 * part / whole


def collect(metacool, mc, layout, stored):
    """The report of the run that wrote the map metacool (see the module's
    docstring), from the map, opened (mc), the layout of the summary of its
    final clusters (layout) and its network (stored), as metator dump reads
    them (see metator.dump.read_run): nothing else is read.

    Returns
    -------
    Report
    """
    metator = mc.metadata.get("metator") if isinstance(mc.metadata.get("metator"), dict) else {}
    entries = {name: metator[name] for name in STEPS if isinstance(metator.get(name), dict)}
    parameters = {name: entry.get("parameters") or {} for name, entry in entries.items()}
    validation = parameters.get("validation", {})
    report = Report()

    # The report itself, and the run.
    report.sections["report"] = {
        "format_version": FORMAT_VERSION,
        "sample": sample_name(metacool),
        "map": str(metacool),
        "metator_version": __version__,
        "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "map_written": mc.info.get("creation-date"),
    }
    inputs = validation.get("input") or validation.get("map")
    partition = parameters.get("partition", {})
    mge = parameters.get("mge", {})
    # The recursive partition of the validation: --rec-iter and --rec-overlap
    # of pipeline, --iterations and --overlap of metator validation.
    recursive = ("rec_iter", "rec_overlap") if "rec_iter" in validation else ("iterations", "overlap")
    run = {
        "steps": ",".join(entries),
        "input": ",".join(map(str, inputs)) if isinstance(inputs, list) else inputs,
        "assembly": validation.get("assembly"),
        "enzyme": parameters.get("network", {}).get("enzyme") or validation.get("enzyme"),
        "depth": parameters.get("network", {}).get("depth") or validation.get("depth"),
        "normalization": (entries.get("partition") or {}).get("normalization") or validation.get("normalization"),
        "algorithm": partition.get("algorithm"),
        "iterations": partition.get("iterations"),
        "overlap_pct": partition.get("overlap"),
        "res_param": partition.get("res_param"),
        "size_bp": validation.get("size"),
        "rec_iter": validation.get(recursive[0]),
        "rec_overlap_pct": validation.get(recursive[1]),
        # --seed of the partition and of the validation (pipeline --seed S:
        # both S), NA for random numbers that no seed set.
        "seed": partition.get("seed"),
        "rec_seed": validation.get("seed"),
    }
    if "mge" in entries:
        run["mge_threshold_bin"] = mge.get("threshold_bin", mge.get("mge_threshold_bin"))
    if mtmc.holds_hosts(mc):
        run["host_min_share_pct"], run["host_min_contigs"] = mth.recorded_thresholds(mc)
    report.sections["run"] = run

    # The contigs, and where they went.
    contigs = mc.contigs
    names = contigs["name"].astype(str).to_numpy()
    lengths = contigs["length"].to_numpy(dtype=np.int64)
    genome_names = mc.genomes["name"].astype(str).to_numpy()
    genome_of = genome_names[contigs["genome_id"].to_numpy()]
    mags = mto.map_mags(mc)
    rounds = [entry.get("name") for entry in mc.cluster_rounds]
    members = {name: mc.cluster_members(name) for name in rounds}
    final_of = dict(zip(members[mtcl.FINAL]["contig"].astype(str), members[mtcl.FINAL]["cluster"].astype(str), strict=True))
    mge_of = (
        dict(zip(members[mtcl.MGE]["contig"].astype(str), members[mtcl.MGE]["cluster"].astype(str), strict=True))
        if mtcl.MGE in members
        else {}
    )
    in_final = np.array([name in final_of for name in names])
    in_mags = np.isin(genome_of, mags)
    in_mges = np.array([name in mge_of for name in names])
    unbinned = genome_of == mtmc.UNBINNED_GENOME
    total = int(lengths.sum())
    n50, l50 = _n50(lengths)
    n90 = int(np.sort(lengths)[::-1][np.argmax(np.cumsum(np.sort(lengths)[::-1]) >= 0.9 * total)]) if total else 0
    gc = pd.to_numeric(contigs["GC_content"], errors="coerce").to_numpy() if "GC_content" in contigs else None
    coverage = pd.to_numeric(contigs["Shotgun_coverage"], errors="coerce").to_numpy() if "Shotgun_coverage" in contigs else None
    gc_pct = None
    if gc is not None and (~np.isnan(gc)).any():
        known = ~np.isnan(gc)
        gc_pct = float(np.average(gc[known], weights=lengths[known]))
    report.sections["contigs"] = {
        "contigs": len(names),
        "contigs_bp": total,
        "contig_n50_bp": n50,
        "contig_l50": l50,
        "contig_n90_bp": n90,
        "contig_max_bp": int(lengths.max()) if len(lengths) else 0,
        "contig_median_bp": float(np.median(lengths)) if len(lengths) else None,
        "gc_pct": gc_pct,
        "coverage_median": float(np.nanmedian(coverage)) if coverage is not None and (~np.isnan(coverage)).any() else None,
        "in_final_clusters_contigs": int(in_final.sum()),
        "in_final_clusters_bp": int(lengths[in_final].sum()),
        "in_final_clusters_pct": _percent(int(lengths[in_final].sum()), total),
        "in_mags_contigs": int(in_mags.sum()),
        "in_mags_bp": int(lengths[in_mags].sum()),
        "in_mags_pct": _percent(int(lengths[in_mags].sum()), total),
        "in_mges_contigs": int(in_mges.sum()),
        "in_mges_bp": int(lengths[in_mges].sum()),
        "unbinned_contigs": int(unbinned.sum()),
        "unbinned_bp": int(lengths[unbinned].sum()),
    }

    # The contacts, and their statistics, as metacooler counted them on the
    # reads (the QC block of the map).
    qc = mc.metadata.get("qc") if isinstance(mc.metadata.get("qc"), dict) else {}
    report.sections["contacts"] = {
        "contacts": qc.get("n_contacts", mc.info.get("sum")),
        "qc_source": qc.get("source"),
        "qc_distance_threshold_bp": qc.get("distance_threshold"),
        "intra_contig_short": qc.get("intra_contig_short"),
        "intra_contig_long": qc.get("intra_contig_long"),
        "inter_contig_intra_genome": qc.get("inter_contig_intra_genome"),
        "inter_genome": qc.get("inter_genome"),
        "trans_contig_ratio": qc.get("trans_contig_ratio"),
        "intra_genome_fraction": qc.get("intra_genome_fraction"),
        "network_contigs": len(stored.table),
        "network_edges": len(stored.network),
        "normalization": stored.normalization,
    }

    # The clusters of each round, and the quality of the final ones, from the
    # values of miComplete (the clusters table holds them at full precision).
    clusters = mc.clusters
    lengths_of = dict(zip(names, lengths, strict=True))
    round_rows = []
    for entry in mc.cluster_rounds:
        name = entry.get("name")
        listed = members[name]
        round_rows.append(
            {
                "round": name,
                "step": entry.get("step"),
                "clusters": int(listed["cluster"].nunique()),
                "contigs": len(listed),
                "length_bp": int(sum(lengths_of.get(contig, 0) for contig in listed["contig"].astype(str))),
            }
        )
    report.tables["rounds"] = pd.DataFrame(round_rows, columns=["round", "step", "clusters", "contigs", "length_bp"])
    final = clusters[(clusters["round"] == mtcl.FINAL).to_numpy()].set_index("name")
    quality = {
        str(name): mtv.bin_quality(row["Weighted completeness"], row["Weighted redundancy"])
        if _number(row.get("Weighted completeness")) is not None and _number(row.get("Weighted redundancy")) is not None
        else None
        for name, row in final.iterrows()
    }
    sizes = final["Length"] if "Length" in final else final["length"]
    by_class = {label: [str(name) for name in final.index if quality[str(name)] == label] for label in QUALITIES}
    counted = {
        "core_clusters": int((clusters["round"] == "core").sum()),
        "overlapping_clusters": int((clusters["round"] == "overlapping").sum()),
        "recursive_steps": sum(1 for name in rounds if str(name).startswith("recursive_")),
        "recursive_clusters": int(clusters["round"].astype(str).str.startswith("recursive_").sum()),
        "final_clusters": len(final),
        "final_clusters_bp": int(pd.to_numeric(sizes, errors="coerce").fillna(0).sum()),
    }
    for label, key in [("HQ", "hq_mags"), ("MQ", "mq_mags"), ("LQ", "lq_mags"), ("Contaminated", "contaminated_mags"),
                       ("Other", "other_clusters")]:  # fmt: skip
        counted[key] = len(by_class[label])
        counted[f"{key}_bp"] = int(pd.to_numeric(sizes[by_class[label]], errors="coerce").fillna(0).sum())
    counted["mags"] = len(mags)
    report.sections["clusters"] = counted

    # The MAGs: the genomes of the map but the MGE genomes, with the values of
    # their validation, and the MGE clusters that they host.
    genomes = mc.genomes.set_index(mc.genomes["name"].astype(str))
    hosted = mto.hosted_mges(mc)
    mag_rows = []
    for name in mags:
        row = genomes.loc[name]
        completeness, redundancy = _number(row.get("Weighted completeness")), _number(row.get("Weighted redundancy"))
        mag_rows.append(
            {
                "name": name,
                "quality": quality.get(name)
                or (mtv.bin_quality(completeness, redundancy) if None not in (completeness, redundancy) else None),
                "length_bp": int(row["length"]),
                "contigs": int(row["n_contigs"]),
                "completeness_pct": None if completeness is None else 100 * completeness,
                "contamination_pct": None if redundancy is None else 100 * (redundancy - 1),
                "gc_pct": _number(row.get("GC-content")),
                "n50_bp": _number(row.get("N50")),
                "hic_abundance": _number(row.get("HiC_abundance")),
                "sg_abundance": _number(row.get("SG_abundance")),
                "intra_genome_ratio": _number(row.get("qc_intra_genome_ratio")),
                "n_mges": len(hosted.get(name, [])),
                "mges": ",".join(hosted.get(name, [])) or None,
            }
        )
    order = {label: position for position, label in enumerate(QUALITIES)}
    mag_rows.sort(key=lambda row: (order.get(row["quality"], len(QUALITIES)), -int(row["length_bp"] or 0), str(row["name"])))
    report.tables["mags"] = pd.DataFrame(mag_rows, columns=list(MAG_COLUMNS))

    # The MGE clusters, and their hosts.
    if mtcl.MGE in members:
        mge_clusters = clusters[(clusters["round"] == mtcl.MGE).to_numpy()].set_index("name")
        mge_genomes = set(mtmc.mge_genomes(mc)) | (set(mge_clusters.index.astype(str)) & set(genome_names))
        losing = {
            final_of[contig]
            for contig, genome in zip(names, genome_of, strict=True)
            if genome in mge_genomes and contig in final_of
        }
        report.sections["mges"] = {
            "mge_clusters": len(mge_clusters),
            "mge_contigs": len(mge_of),
            "mge_bp": int(sum(lengths_of.get(contig, 0) for contig in mge_of)),
            "mge_genomes": len(mge_genomes),
            "mags_losing_contigs": len(losing),
            "mge_random_clusters": int((clusters["round"] == mtcl.MGE_RANDOM).sum()) if mtcl.MGE_RANDOM in members else None,
        }
        mge_rows, links = [], []
        for name, row in mge_clusters.iterrows():
            listed = [contig for contig, cluster in mge_of.items() if cluster == str(name)]
            hosts = _names(row.get("hosts"))
            shares = [float(share) for share in _names(row.get("host_shares"))]
            links += [{"mge": str(name), "mag": host, "share_pct": share} for host, share in zip(hosts, shares, strict=False)]
            mge_rows.append(
                {
                    "name": str(name),
                    "contigs": len(listed),
                    "length_bp": int(sum(lengths_of.get(contig, 0) for contig in listed)),
                    "gc_pct": _number(row.get("GC")),
                    "n_hosts": len(hosts),
                    "hosts": ",".join(hosts) or None,
                    "host_shares_pct": ",".join(f"{share:.2f}" for share in shares) or None,
                }
            )
        report.tables["mge_clusters"] = pd.DataFrame(
            mge_rows, columns=["name", "contigs", "length_bp", "gc_pct", "n_hosts", "hosts", "host_shares_pct"]
        )
        if mtmc.holds_hosts(mc):
            n_hosts = np.array([row["n_hosts"] for row in mge_rows], dtype=int)
            report.sections["hosts"] = {
                "min_share_pct": run.get("host_min_share_pct"),
                "min_contigs": run.get("host_min_contigs"),
                "mge_clusters_with_host": int((n_hosts > 0).sum()),
                "mge_clusters_with_one_host": int((n_hosts == 1).sum()),
                "mge_clusters_with_several_hosts": int((n_hosts > 1).sum()),
                "mge_clusters_without_host": int((n_hosts == 0).sum()),
                "mags_hosting_mges": sum(1 for name in mags if hosted.get(name)),
            }
            report.tables["hosts"] = pd.DataFrame(links, columns=["mge", "mag", "share_pct"])
        else:
            # No host run yet: the MGE clusters have no host to count.
            report.tables["mge_clusters"]["n_hosts"] = None

    # The command lines of the steps.
    report.sections["commands"] = {name: entry.get("command") for name, entry in entries.items()}

    # What the plots of report.html are drawn from.
    report.plots = _plot_data(report, names, lengths, genome_of, mags, mge_of, gc, coverage, contigs, stored)
    return report


def _plot_data(report, names, lengths, genome_of, mags, mge_of, gc, coverage, contigs, stored):
    """What the plots of report.html are drawn from (see collect): the contigs
    (at most MAX_PLOTTED_CONTIGS), their GC content and coverage (the shotgun
    coverage, or the Hi-C hits per kb without it), and the contacts between
    the MAGs, summed over the network of the map (at most MAX_HEATMAP_MAGS)."""
    category = np.where(np.isin(genome_of, mags), "MAG", np.where(genome_of == mtmc.UNBINNED_GENOME, "unbinned", "MGE"))
    category = np.where([name in mge_of for name in names], "MGE", category)
    if coverage is not None and (~np.isnan(coverage)).any():
        depth, depth_label = coverage, "Shotgun coverage"
    else:
        hits = pd.to_numeric(contigs["Hit"], errors="coerce").to_numpy() if "Hit" in contigs else np.full(len(names), np.nan)
        depth, depth_label = 1000 * hits / np.maximum(lengths, 1), "Hi-C hits per kb"
    shown = np.argsort(np.where(category == "unbinned", 1, 0) * 10**15 - lengths, kind="stable")[:MAX_PLOTTED_CONTIGS]
    plotted = {
        "name": names[shown].tolist(),
        "length": lengths[shown].tolist(),
        "genome": genome_of[shown].tolist(),
        "category": category[shown].tolist(),
        "gc": None if gc is None else [None if math.isnan(value) else float(value) for value in gc[shown]],
        "depth": [None if math.isnan(value) else float(value) for value in np.asarray(depth, dtype=float)[shown]],
        "depth_label": depth_label,
        "shown": len(shown),
        "total": len(names),
    }
    # The distribution of the lengths of the contigs, by category: their
    # number in bins of 0.1 log10 bp, binned here, so that the page holds a
    # count per bin, not a length per contig.
    logs = np.log10(np.maximum(lengths, 1))
    low, high = (int(np.floor(logs.min() * 10)), int(np.ceil(logs.max() * 10))) if len(logs) else (0, 1)
    edges = np.arange(low, max(high, low + 1) + 1) / 10
    distribution = {
        "edges": edges.tolist(),
        "counts": {label: np.histogram(logs[category == label], bins=edges)[0].tolist() for label in ("MAG", "MGE", "unbinned")},
    }
    # The contacts between the MAGs: the weights of the network, summed by
    # MAG, the largest MAGs of each class first.
    table = report.tables["mags"]
    heatmap_mags = table["name"].tolist()[:MAX_HEATMAP_MAGS]
    position = {name: index for index, name in enumerate(heatmap_mags)}
    mag_of_contig = np.full(len(names) + 1, -1)
    for index, genome in enumerate(genome_of):
        mag_of_contig[index + 1] = position.get(genome, -1)
    network = stored.network
    first = mag_of_contig[network["contig1"].to_numpy()]
    second = mag_of_contig[network["contig2"].to_numpy()]
    kept = (first >= 0) & (second >= 0)
    matrix = np.zeros((len(heatmap_mags), len(heatmap_mags)))
    np.add.at(matrix, (first[kept], second[kept]), network["weight"].to_numpy()[kept])
    matrix = matrix + matrix.T - np.diag(np.diag(matrix))
    return {
        "contigs": plotted,
        "lengths": distribution,
        "heatmap": {"mags": heatmap_mags, "contacts": matrix.tolist(), "normalization": stored.normalization},
    }


# ---------------------------------------------------------------------------
# report.txt
# ---------------------------------------------------------------------------


def _text(key, value):
    """A value of report.txt as text: NA for a missing one, an int as an int, a
    float with 2 decimals for a percent (_pct), 4 for a ratio or a fraction,
    else as Python writes it shortest; a text on one line."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "true" if value else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if float(value).is_integer() and not key.endswith(("_pct", "_ratio", "_fraction")):
            return str(int(value))
        if key.endswith("_pct"):
            return f"{value:.2f}"
        if key.endswith(("_ratio", "_fraction")):
            return f"{value:.4f}"
        return repr(float(value))
    return " ".join(str(value).split()) or "NA"


def write_txt(report, path):
    """Write report.txt at path (see the module's docstring): a header of
    comments; each section, ``[name]``, then a line per key, the key, spaces,
    and its value; each table, ``[table: name]``, then its header and a line
    per row, their fields aligned on spaces. No key, column or field holds a
    space, but the values of the sections run and commands (the rest of their
    line): a line splits on its first spaces. A missing value is NA."""
    lines = [
        f"# metaTOR report of {report.sections['report']['sample']}",
        f"# format {FORMAT_VERSION}: [section] then a key and its value per line; [table: name] then a header and a row",
        "# per line, on spaces; NA: missing. See https://koszullab.github.io/metator/formats/#the-report-of-a-run",
        "",
    ]
    for name, values in report.sections.items():
        width = max((len(key) for key in values), default=0)
        lines.append(f"[{name}]")
        lines += [f"{key.ljust(width)}  {_text(key, value)}" for key, value in values.items()]
        lines.append("")
    for name, table in report.tables.items():
        lines.append(f"[table: {name}]")
        rows = [[_text(column, value).replace(" ", "_") for column, value in zip(table.columns, row, strict=True)]
                for row in table.itertuples(index=False)]  # fmt: skip
        widths = [max([len(column)] + [len(row[index]) for row in rows]) for index, column in enumerate(table.columns)]
        lines.append("  ".join(column.ljust(width) for column, width in zip(table.columns, widths, strict=True)).rstrip())
        lines += ["  ".join(field.ljust(width) for field, width in zip(row, widths, strict=True)).rstrip() for row in rows]
        lines.append("")
    with open(path, "w") as handle:
        handle.write("\n".join(lines))


def read_txt(path):
    """Read report.txt (see write_txt) back: its sections, by name, each its
    keys and their values, as text; and its tables, by name, as
    pandas.DataFrame of text, NA missing. What a MultiQC module or a script
    does to aggregate the reports of runs."""
    sections, tables, current, header, rows = {}, {}, None, None, []

    def close():
        if header is not None:
            tables[current] = pd.DataFrame(rows, columns=header, dtype=object)

    with open(path) as handle:
        for line in handle.read().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            match = re.fullmatch(r"\[(table: )?([a-z_]+)\]", line.strip())
            if match:
                close()
                current, header, rows = match.group(2), None, []
                if match.group(1):
                    header = []
                else:
                    sections[current] = {}
                continue
            if header == []:
                header = line.split()
            elif header is not None:
                rows.append([None if field == "NA" else field for field in line.split()])
            else:
                key, _, value = line.partition(" ")
                value = value.strip()
                sections[current][key] = None if value == "NA" else value
    close()
    return sections, tables


# ---------------------------------------------------------------------------
# report_mqc.yaml
# ---------------------------------------------------------------------------


def write_mqc(report, path):
    """Write report_mqc.yaml at path: MultiQC's custom content, a row of its
    General Statistics table for the run (MQC_COLUMNS), named after it, which
    MultiQC aggregates over the runs whose reports it finds (``multiqc .``):
    the same id for every run."""
    import yaml

    sample = report.sections["report"]["sample"]
    row = {}
    pconfig = []
    for key, section, config in MQC_COLUMNS:
        value = report.sections.get(section, {}).get(key)
        if value is None:
            continue
        row[key] = int(value) if isinstance(value, (int, np.integer)) else float(value)
        pconfig.append({key: dict(config)})
    content = {
        "id": "metator",
        "section_name": "metaTOR",
        "description": f"metaTOR runs: the MAGs of the map of each run, and its MGE clusters (format {FORMAT_VERSION}).",
        "plot_type": "generalstats",
        "pconfig": pconfig,
        "data": {sample: row},
    }
    with open(path, "w") as handle:
        handle.write(f"# metaTOR report of {sample}: MultiQC custom content (General Statistics); multiqc finds it.\n")
        yaml.safe_dump(content, handle, sort_keys=False, default_flow_style=False, allow_unicode=True)


# ---------------------------------------------------------------------------
# report.html
# ---------------------------------------------------------------------------


def _json_ready(value):
    """value with NaN and infinities as None, and numpy's types as Python's: what
    JSON.parse reads (it refuses NaN)."""
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def figures(report):
    """The plots of report.html, by the id of their place in the page: Plotly
    figures (data and layout), which plotly.js draws."""
    sections, tables, plots = report.sections, report.tables, report.plots
    figs = {}
    font = {"family": "Inter, Segoe UI, Helvetica, Arial, sans-serif", "size": 12}
    # A log axis: its labels at the powers of 10, a gridline at each digit.
    log_axis = {"type": "log", "dtick": 1, "minor": {"dtick": "D1", "showgrid": True, "gridcolor": "#f0f0f0"}}
    # The ticks of the axis of the lengths of the contigs (log10 bp): 1, 2, 5.
    length_ticks = [(100, "100"), (200, "200"), (500, "500"), (1e3, "1 k"), (2e3, "2 k"), (5e3, "5 k"), (1e4, "10 k"),
                    (2e4, "20 k"), (5e4, "50 k"), (1e5, "100 k"), (2e5, "200 k"), (5e5, "500 k"), (1e6, "1 M"),
                    (2e6, "2 M"), (5e6, "5 M"), (1e7, "10 M")]  # fmt: skip
    base = {"font": font, "margin": {"l": 60, "r": 20, "t": 30, "b": 50}, "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(0,0,0,0)", "hovermode": "closest"}  # fmt: skip

    # The quality of the MAGs: completeness against contamination, a dot per
    # MAG, its size its length.
    mags = tables["mags"]
    data = []
    for label in QUALITIES:
        part = mags[mags["quality"] == label]
        if not len(part):
            continue
        data.append(
            {
                "type": "scatter",
                "mode": "markers",
                "name": f"{label} ({len(part)})",
                "x": part["completeness_pct"].tolist(),
                "y": part["contamination_pct"].tolist(),
                "text": part["name"].tolist(),
                "customdata": list(
                    zip(part["length_bp"].tolist(), part["contigs"].tolist(), part["n_mges"].tolist(), strict=True)
                ),
                "marker": {
                    "color": QUALITY_COLOURS[label],
                    "size": (6 + 3 * np.sqrt(part["length_bp"].to_numpy() / 1e5)).clip(6, 40).tolist(),
                    "line": {"width": 0.5, "color": "#333"},
                    "opacity": 0.85,
                },
                "hovertemplate": "<b>%{text}</b><br>completeness %{x:.1f}%<br>contamination %{y:.1f}%<br>"
                "%{customdata[0]:,} bp, %{customdata[1]} contigs<br>%{customdata[2]} MGE cluster(s)<extra></extra>",
            }
        )
    shapes = [
        {"type": "line", "x0": threshold, "x1": threshold, "y0": 0, "y1": 1, "yref": "paper", "line": {"dash": "dot", "color": "#999", "width": 1}}
        for threshold in (50, 70, 90)
    ] + [
        {"type": "line", "x0": 0, "x1": 1, "xref": "paper", "y0": threshold, "y1": threshold, "line": {"dash": "dot", "color": "#999", "width": 1}}
        for threshold in (5, 10)
    ]  # fmt: skip
    figs["plot-quality"] = {
        "data": data,
        "layout": {**base, "xaxis": {"title": "Completeness (%)", "range": [-2, 105]},
                   "yaxis": {"title": "Contamination (%)",
                             "range": [-1, max(12.0, 1.1 * float(np.nanmax(mags["contamination_pct"].to_numpy(dtype=float, na_value=np.nan), initial=0)))]},
                   "shapes": shapes, "height": 460},
    }  # fmt: skip

    # The final clusters of each class: how many, and their length.
    clusters = sections["clusters"]
    keys = ["hq_mags", "mq_mags", "lq_mags", "contaminated_mags", "other_clusters"]
    lengths = [clusters[f"{key}_bp"] / 1e6 for key in keys]
    figs["plot-classes"] = {
        "data": [
            {"type": "bar", "x": list(QUALITIES), "y": [clusters[key] for key in keys], "customdata": lengths,
             "text": [f"{clusters[key]} · {length:.2f} Mbp" for key, length in zip(keys, lengths, strict=True)],
             "textposition": "outside", "cliponaxis": False,
             "marker": {"color": [QUALITY_COLOURS[label] for label in QUALITIES]},
             "hovertemplate": "%{x}: %{y} cluster(s), %{customdata:.2f} Mbp<extra></extra>"},
        ],
        "layout": {**base, "height": 320, "yaxis": {"title": "Final clusters", "rangemode": "tozero"}},
    }  # fmt: skip

    # The clusters of each round of the clustering.
    rounds = tables["rounds"]
    figs["plot-rounds"] = {
        "data": [{"type": "bar", "x": rounds["round"].tolist(), "y": rounds["clusters"].tolist(),
                  "text": rounds["clusters"].tolist(), "textposition": "outside", "cliponaxis": False, "marker": {"color": "#4c78a8"},
                  "customdata": rounds["length_bp"].tolist(),
                  "hovertemplate": "%{x}: %{y} clusters, %{customdata:,} bp<extra></extra>"}],
        "layout": {**base, "height": 320, "yaxis": {"title": "Clusters", "rangemode": "tozero"}},
    }  # fmt: skip

    # The length of the contigs: in a MAG, in an MGE cluster, in none, by
    # bins of 0.1 log10 bp (see _plot_data).
    colours = {"MAG": "#4c78a8", "MGE": "#f58518", "unbinned": "#bab0ac"}
    edges = plots["lengths"]["edges"]
    centres = [(left + right) / 2 for left, right in itertools.pairwise(edges)]
    ranges = [[10**left, 10**right] for left, right in itertools.pairwise(edges)]
    figs["plot-lengths"] = {
        "data": [{"type": "bar", "name": label, "x": centres, "y": [count or None for count in counts], "width": 0.1,
                  "opacity": 0.75, "marker": {"color": colours[label]}, "customdata": ranges,
                  "hovertemplate": f"{label}: %{{y}} contigs of %{{customdata[0]:,.0f}} to %{{customdata[1]:,.0f}} bp<extra></extra>"}
                 for label, counts in plots["lengths"]["counts"].items() if any(counts)],
        "layout": {**base, "barmode": "overlay", "height": 320,
                   "xaxis": {"title": "Contig length (bp)", "tickvals": [math.log10(value) for value, _ in length_ticks],
                             "ticktext": [label for _, label in length_ticks]},
                   "yaxis": {"title": "Contigs", **log_axis}},
    }  # fmt: skip

    # The contigs, their GC content against their coverage, coloured by MAG.
    contigs = plots["contigs"]
    if contigs["gc"] is not None:
        frame = pd.DataFrame({key: contigs[key] for key in ("name", "length", "genome", "category", "gc", "depth")})
        largest = mags["name"].tolist()[:20]
        groups = [(name, frame[frame["genome"] == name]) for name in largest]
        groups.append(("other MAGs", frame[(frame["category"] == "MAG") & ~frame["genome"].isin(largest)]))
        groups.append(("MGE clusters", frame[frame["category"] == "MGE"]))
        groups.append(("unbinned", frame[frame["category"] == "unbinned"]))
        data = []
        for name, part in groups:
            if not len(part):
                continue
            trace = {
                "type": "scattergl", "mode": "markers", "name": name, "x": part["gc"].tolist(), "y": part["depth"].tolist(),
                "text": part["name"].tolist(), "customdata": list(zip(part["length"].tolist(), part["genome"].tolist(), strict=True)),
                "marker": {"size": np.clip(2 + np.log10(np.maximum(part["length"].to_numpy(), 1)) - 2, 2, 9).tolist(), "opacity": 0.7},
                "hovertemplate": "<b>%{text}</b><br>%{customdata[1]}<br>%{customdata[0]:,} bp<br>GC %{x:.1f}%<br>"
                + contigs["depth_label"] + " %{y:.2f}<extra></extra>",
            }  # fmt: skip
            if name == "unbinned":
                trace["marker"]["color"] = "#cccccc"
            elif name == "MGE clusters":
                trace["marker"]["color"] = "#f58518"
            elif name == "other MAGs":
                trace["marker"]["color"] = "#7f7f7f"
            data.append(trace)
        shown = "" if contigs["shown"] == contigs["total"] else f" (the {contigs['shown']} longest of {contigs['total']})"
        figs["plot-contigs"] = {
            "data": data,
            "layout": {**base, "height": 520, "xaxis": {"title": "GC content (%)"},
                       "yaxis": {"title": contigs["depth_label"], **log_axis},
                       "title": {"text": f"Contigs{shown}", "font": {"size": 13}}},
        }  # fmt: skip

    # The contacts of the reads: within a contig, short and long range, between
    # the contigs of a genome, between genomes.
    contacts = sections["contacts"]
    categories = [("intra_contig_short", "Within a contig, < 1 kb"), ("intra_contig_long", "Within a contig, ≥ 1 kb"),
                  ("inter_contig_intra_genome", "Between contigs of a genome"), ("inter_genome", "Between genomes")]  # fmt: skip
    if all(contacts.get(key) is not None for key, _ in categories):
        figs["plot-contacts"] = {
            "data": [{"type": "bar", "orientation": "h", "y": [label for _, label in categories],
                      "x": [contacts[key] for key, _ in categories], "marker": {"color": ["#9ecae1", "#4292c6", "#1a9850", "#d73027"]},
                      "hovertemplate": "%{y}: %{x:,} contacts<extra></extra>"}],
            "layout": {**base, "height": 260, "margin": {**base["margin"], "l": 200}, "xaxis": {"title": "Contacts"}},
        }  # fmt: skip

    # The contacts between the MAGs (the network of the map, summed).
    heatmap = plots["heatmap"]
    if len(heatmap["mags"]) > 1:
        values = np.asarray(heatmap["contacts"], dtype=float)
        logged = np.where(values > 0, np.log10(np.where(values > 0, values, 1)), np.nan)
        figs["plot-heatmap"] = {
            "data": [{"type": "heatmap", "z": logged.tolist(), "x": heatmap["mags"], "y": heatmap["mags"],
                      "colorscale": "Reds", "colorbar": {"title": "log10"},
                      "hovertemplate": "%{y} with %{x}<br>log10 weight %{z:.2f}<extra></extra>"}],
            "layout": {**base, "height": 640, "margin": {"l": 170 if len(heatmap["mags"]) <= 40 else 40, "r": 20, "t": 40,
                                                         "b": 150 if len(heatmap["mags"]) <= 40 else 40},
                       "xaxis": {"showticklabels": len(heatmap["mags"]) <= 40, "tickangle": -45},
                       "yaxis": {"showticklabels": len(heatmap["mags"]) <= 40, "autorange": "reversed"},
                       "title": {"text": f"{heatmap['normalization']} weights of the network, summed by MAG", "font": {"size": 12}}},
        }  # fmt: skip

    # The hosts of the MGE clusters: from each MGE cluster to its host MAGs.
    if "hosts" in tables and len(tables["hosts"]):
        links = tables["hosts"]
        nodes = list(dict.fromkeys(links["mge"].tolist())) + list(dict.fromkeys(links["mag"].tolist()))
        index = {node: position for position, node in enumerate(nodes)}
        n_mges = links["mge"].nunique()
        figs["plot-hosts"] = {
            "data": [{"type": "sankey", "arrangement": "snap",
                      "node": {"label": nodes, "pad": 8, "thickness": 14,
                               "color": ["#f58518"] * n_mges + ["#4c78a8"] * (len(nodes) - n_mges)},
                      "link": {"source": [index[node] for node in links["mge"]], "target": [index[node] for node in links["mag"]],
                               "value": links["share_pct"].tolist(),
                               # Plotly draws them white on a transparent paper.
                               "color": "rgba(245, 133, 24, 0.35)",
                               "hovertemplate": "%{source.label} → %{target.label}: %{value:.1f}% of its signal<extra></extra>"}}],
            "layout": {**base, "height": max(320, 22 * len(nodes))},
        }  # fmt: skip
    return figs


def _links(outdir, report):
    """The outputs of the run that report.html links to, which are in its output
    directory, outdir, as metator export writes them: by MAG, its FASTA files
    and its cool file; the tables and the plots."""

    def found(*parts):
        path = join(*parts)
        return path if os.path.exists(join(outdir, path)) else None

    mags = {
        name: {"fasta": found(mto.MAGS_DIR, f"{name}.fa"), "stitched": found(mto.MAGS_DIR, f"{name}{mto.STITCHED}"),
               "cool": found(mto.COOL_DIR, f"{name}.cool")}
        for name in report.tables["mags"]["name"]
    }  # fmt: skip
    files = [path for path in (
        [join(mto.TABLES_DIR, name) for name in sorted(os.listdir(join(outdir, mto.TABLES_DIR)))] if os.path.isdir(join(outdir, mto.TABLES_DIR)) else []
    ) + [join(mtf.PLOT_DIR, name) for name in (*mtf.PLOTS, mth.PLOT_FILE)] + [REPORT_TXT, REPORT_MQC] if os.path.exists(join(outdir, path))]  # fmt: skip
    return {"mags": mags, "files": files}


def write_html(report, path, outdir):
    """Write report.html at path, from its Jinja2 template: the run, its numbers
    and tables, and its plots (see figures), which plotly.js draws once the
    page has loaded it from its CDN (PLOTLY_JS); the outputs of the run that
    are in its output directory, outdir, as links."""
    import jinja2

    environment = jinja2.Environment(
        loader=jinja2.PackageLoader("metator", "templates"),
        autoescape=jinja2.select_autoescape(["html", "j2"]),
        trim_blocks=True,
        lstrip_blocks=True,
        undefined=jinja2.StrictUndefined,
    )
    environment.filters["value"] = lambda value, key="": _text(key, value)
    template = environment.get_template("report.html.j2")
    tables = {name: _json_ready(table.to_dict(orient="records")) for name, table in report.tables.items()}
    page = template.render(
        report=report.sections,
        tables=tables,
        columns={name: list(table.columns) for name, table in report.tables.items()},
        figures=_json_ready(figures(report)),
        links=_links(outdir, report),
        qualities=QUALITIES,
        colours=QUALITY_COLOURS,
        plotly_js=PLOTLY_JS,
        format_version=FORMAT_VERSION,
    )
    with open(path, "w") as handle:
        handle.write(page)


def write_report(report, outdir):
    """Write the report of a run (see collect) in outdir (see the module's
    docstring): report.txt, report_mqc.yaml and report.html, which replace
    files of their names. The map is not read again: what report.html links
    to is what outdir holds as it is written.

    Returns
    -------
    list of str
        Their paths, in the order of REPORT_FILES."""
    os.makedirs(outdir, exist_ok=True)
    paths = [join(outdir, name) for name in REPORT_FILES]
    for path in paths:
        if os.path.islink(path):
            os.remove(path)
    # report.html last, which links to the other two.
    write_txt(report, paths[0])
    write_mqc(report, paths[2])
    write_html(report, paths[1], outdir)
    return paths
