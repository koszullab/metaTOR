import sys
import os
import metator.io as mio
import metator.host as mth
import networkx as nx
import numpy as np
import pandas as pd
import pytest # type: ignore
import tempfile

#LISTE DES PARAMETRES D'ENTRÉES DU MODULE

contig_data_file = "tests_data/host_test_data/contig_data_final.txt"
binning_file = "tests_data/host_test_data/bin_summary.txt"
network_file =  "tests_data/host_test_data/network_0.txt"
mges_bin_summary_file = "tests_data/host_test_data/mges_bin_summary.tsv"

threshold_association = 10 #(default)
interacting_contig = 5  #(default)

contig_data = pd.read_csv(contig_data_file, sep="\t")
bin_summary = pd.read_csv(binning_file, sep="\t")
network_data = pd.read_csv(network_file, sep="\t", names=["contig1", "contig2", "signal"])
mges_bin_summary = pd.read_csv(mges_bin_summary_file, sep="\t")



def test_create_bins():
    mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
    assert isinstance(mags, dict)
    assert isinstance(mge_mags, dict)
    assert len(mags) == 50
    assert len(mge_mags) == 566

    # Vérifications spécifiques
    assert "metator_00010_00000" in mags
    assert "MetaTOR_MGE_00010" in mge_mags

    mag = mags["metator_00010_00000"]
    mge = mge_mags["MetaTOR_MGE_00010"]

    assert isinstance(mag, mth.Mag)
    assert isinstance(mge, mth.MgeMag)

    assert mag.contigs["NODE_5489_len_50101"] == 152
    assert mag.contigs["NODE_19150_len_13073"] == 320
    assert mge.contigs["NODE_117_len_4104"] == 721
    assert mge.contigs["NODE_6974_len_5718"] == 10117
    

def test_map_contigs_to_mags():
    mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
    contig_to_mag = mth.map_contigs_to_mags(mags)
    assert isinstance(contig_to_mag, dict)
    assert len(contig_to_mag) == 7006

    # # Vérifie le mapping correct
    assert contig_to_mag[10435] == "metator_00005_00004"
    assert contig_to_mag[10669] == "metator_00005_00001"


def test_build_contig_graph():
    G = mth.build_contig_graph(network_data)
    assert isinstance(G, nx.Graph) 

def test_estimate_noise():
    
    mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
    contig_to_mag = mth.map_contigs_to_mags(mags)
    # Appel de la fonction
    result = mth.estimate_noise(mags, network_data)

    # Vérifications
    assert isinstance(mags["metator_00005_00004"].intra_signal, float)
    assert mags["metator_00005_00004"].intra_signal >= 0
    assert len(mags["metator_00005_00004"].inter_signals) == len(mags) - 1

    assert isinstance(mags["metator_00005_00001"].intra_signal, float)
    assert mags["metator_00005_00001"].intra_signal >= 0
    assert len(mags["metator_00005_00001"].inter_signals) == len(mags) - 1

    assert result == mags





# def test_plot_mag_noise():
#     mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
#     mags = mth.estimate_noise(mags, network_data)
#     #mags = mth.classify_mags(mags, bin_summary)

#     # Injecte une interaction simulée pour que le plot fonctionne
#     # for mag in mags.values():
#     #     mag.inter_signals = {"metator_7000_00000": 1.23}  # fake MAG interaction

#     with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_file:
#         image_path = tmp_file.name

#     try:
#         mth.plot_mag_noise(mags, image_file=image_path)
#         assert os.path.exists(image_path)
#         assert os.path.getsize(image_path) > 0
#         print("✅ test_plot_mag_noise passed.")
#     finally:
#         if os.path.exists(image_path):
#             os.remove(image_path)


def test_classify_mags():
    mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
    bin_summary.rename(columns={"Unnamed: 0": "MAG"}, inplace=True)
    mag_class = mth.classify_mags(mags, bin_summary)
    assert isinstance(mag_class["metator_00005_00001"].quality, str)
    assert mags["metator_00005_00001"].quality == "MQ"
   

def test_compute_mge_mag_interactions():

    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = os.path.join(tmpdir, "test_output.tsv")
        image_file = os.path.join(tmpdir, "test_plot.png")

        mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)

        updated_mags, updated_mge_mags = mth.compute_mge_mag_interactions(
            network_data, mags, mge_mags,
            interaction_threshold=10.0,
            min_interacting_contigs=5,
            output_file=output_file,
            image_file=image_file
        )

        # Vérifie que le fichier a été créé
        assert os.path.isfile(output_file)
        assert os.path.isfile(image_file)

        # Vérifie que les interactions ont été bien enregistrées
        with open(output_file, "r") as f:
            lines = f.readlines()
            assert len(lines) >= 2  # au moins un header + une interaction

        # Vérifie que la liaison a été faite entre MGE1 et MAG1/MAG2
        assert "metator_00061_00000" in updated_mge_mags["MetaTOR_MGE_00042"].hosts or "metator_00013_00000" in updated_mge_mags["MetaTOR_MGE_00103"].hosts
        assert "MetaTOR_MGE_00042" in updated_mags["metator_00061_00000"].mges or "MetaTOR_MGE_00185" in updated_mags["metator_00061_00000"].mges


def test_annotate_hosts():

    bin_summary.rename(columns={"Unnamed: 0": "MAG"}, inplace=True)

     # Run the function
    mags, mge_mags = mth.annotate_hosts(
        contig_data,
        network_data,
        bin_summary,
        mges_bin_summary,
        interaction_threshold=10,
        min_interacting_contigs=5
    )

    # Check that the return values are as expected
    assert isinstance(mags, dict)
    assert isinstance(mge_mags, dict)
    assert "metator_00061_00000" in mags
    assert "MetaTOR_MGE_00196" in mge_mags
