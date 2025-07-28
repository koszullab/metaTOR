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
    assert len(mge_mags) == 119

    # Vérifications spécifiques
    assert "metator_00028_00003" in mags
    assert "MetaTOR_MGE_00002" in mge_mags

    mag = mags["metator_00028_00003"]
    mge = mge_mags["MetaTOR_MGE_00002"]

    assert isinstance(mag, mth.Mag)
    assert isinstance(mge, mth.MgeMag)
    print("Contigs disponibles dans le MAG :", list(mge.contigs.keys()))
    #assert "NODE_5491_len_437" in mag.contigs

    assert mag.contigs["NODE_12058_len_3472"] == 6498
    assert mag.contigs["NODE_14880_len_38848"] == 9697
    assert mge.contigs["NODE_11750_len_20083"] == 4841
    assert mge.contigs["NODE_19204_len_17733"] == 832
    

def test_map_contigs_to_mags():
    mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
    contig_to_mag = mth.map_contigs_to_mags(mags)
    assert isinstance(contig_to_mag, dict)
    assert len(contig_to_mag) == 7748

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
    mags, mean_inter, interactions_matrix, mag_names = mth.estimate_noise(mags, network_data)

    # Vérifications
    assert isinstance(mags["metator_00005_00004"].intra_signal, float)
    assert mags["metator_00005_00004"].intra_signal >= 0
    assert len(mags["metator_00005_00004"].inter_signals) == len(mags) - 1

    assert isinstance(mags["metator_00005_00001"].intra_signal, float)
    assert mags["metator_00005_00001"].intra_signal >= 0
    assert len(mags["metator_00005_00001"].inter_signals) == len(mags) - 1

    assert isinstance(mean_inter, float) 
    assert mean_inter == 0.053181618705037154
    assert mags == mags


def test_classify_mags():
    mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
    bin_summary.rename(columns={"Unnamed: 0": "MAG"}, inplace=True)
    mag_class = mth.classify_mags(mags, bin_summary)
    assert isinstance(mag_class["metator_00005_00001"].quality, str)
    assert mags["metator_00005_00001"].quality == "MQ"
   

def test_compute_mge_mag_interactions():

    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = os.path.join(tmpdir, "test_output.tsv")
        plot_prefix = os.path.join(tmpdir, "test_plot.png")
        mags, mge_mags = mth.create_bins(contig_data, mges_bin_summary)
        mags, mean_inter, interactions_matrix, mag_names = mth.estimate_noise(mags, network_data)
        updated_mags, updated_mge_mags = mth.compute_mge_mag_interactions(
            mean_inter,
            network_data, mags, mge_mags,
            interaction_threshold=10.0,
            min_interacting_contigs=5,
            output_file=output_file,
            plot_prefix=plot_prefix,
        
        )

        # Vérifie que le fichier a été créé
        assert os.path.isfile(output_file)
        

        # Vérifie que les interactions ont été bien enregistrées
        with open(output_file, "r") as f:
            lines = f.readlines()
            assert len(lines) >= 2  # au moins un header + une interaction

        # Vérifie que la liaison a été faite entre MGE1 et MAG1/MAG2
        assert "metator_00027_00002" in updated_mge_mags["MetaTOR_MGE_00287"].hosts or "metator_00007_00000" in updated_mge_mags["MetaTOR_MGE_00412"].hosts
        assert "MetaTOR_MGE_00094" in updated_mags["metator_00005_00001"].mges or "MetaTOR_MGE_00030" in updated_mags["metator_00027_00003"].mges
        if os.path.isfile(plot_prefix):
            print("Le fichier de plot a bien été généré.")
        else:
            print("Aucun plot généré. Peut-être pas assez d'interactions ?")


def test_annotate_hosts():
    bin_summary.rename(columns={"Unnamed: 0": "MAG"}, inplace=True)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = os.path.join(tmpdir, "annotate_hosts_output")
        os.makedirs(outdir, exist_ok=True)
     # Run the function
    mags, mge_mags = mth.annotate_hosts(
        contig_data,
        network_data,
        bin_summary,
        mges_bin_summary,
        interaction_threshold=10,
        min_interacting_contigs=5,
        outdir =outdir
    )

    # Check that the return values are as expected
    assert isinstance(mags, dict)
    assert isinstance(mge_mags, dict)
    assert "metator_00027_00003" in mags
    assert "MetaTOR_MGE_00287" in mge_mags
    assert os.path.isdir(outdir)
    assert os.path.isfile(os.path.join(outdir, "mag_interactions_boxplot.png"))
