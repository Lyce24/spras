import warnings
from pathlib import Path

import pandas as pd

from src.prm import PRM
from src.util import add_rank_column, prepare_volume, run_container

__all__ = ["TieDIE"]

class TieDIE(PRM):
    # we need edges (weighted), source set (with prizes), and target set (with prizes).
    required_inputs = ["edges", "sources", "targets"]

    @staticmethod
    def generate_inputs(data, filename_map):
        """
        Access fields from the dataset and write the required input files
        @param data: dataset
        @param filename_map: a dict mapping file types in the required_inputs to the filename for that type
        """
        # ensures the required input are within the filename_map
        for input_type in TieDIE.required_inputs:
            if input_type not in filename_map:
                raise ValueError(f"{input_type} filename is missing")

        # will take the sources and write them to files, and repeats with targets
        for node_type in ["sources", "targets"]:
            nodes = data.request_node_columns([node_type])
            # check if the nodes have prizes or not
            if data.contains_node_columns("prize"):
                node_df = data.request_node_columns(["prize"])
                nodes = pd.merge(nodes, node_df, on="NODEID")
                nodes["sign"] = "+"
                # creates with the node type without headers
                nodes.to_csv(filename_map[node_type],index=False,sep="\t",columns=["NODEID", "prize", "sign"],header=False)
            else:
                # If there aren't prizes but are sources and targets, make prizes based on them
                nodes = data.request_node_columns([node_type])
                # make all nodes have a prize of 1
                nodes["prize"] = 1.0
                nodes["sign"] = "+"
                # creates with the node type without headers
                nodes.to_csv(filename_map[node_type],index=False,sep="\t",columns=["NODEID", "prize", "sign"],header=False)

        # create the network of edges
        edges = data.get_interactome()
        edges["type"] = "-a>"
        # drop the weight column
        edges = edges.drop(columns=["Weight"])
        # creates the edges files that contains the head and tail nodes and the weights after them
        edges.to_csv(filename_map["edges"],sep="\t",index=False,columns=["Interactor1", "type", "Interactor2"],header=False)

    # Skips parameter validation step
    @staticmethod
    def run(edges=None,sources=None,targets=None,output_file=None,s : float = 1.0,c : int = 3, p : int = 1000, pagerank : bool = False, all_paths : bool = False, singularity=False):
        """
        Run TieDIE with Docker
        @param source:  input node types with sources (required)
        @param target:  input node types with targets (required)
        @param network:  input network file (required)
        @param output_file: path to the output pathway file (required)
        @param s: Network size control factor (optional) (default 1)
        @param d_expr: List of significantly differentially expressed genes, along with log-FC or FC values (i.e. by edgeR for RNA-Seq or SAM for microarray data. Generated by a sample-dichotomy of interest. (optional)
        @param a: Linker Cutoff (overrides the Size factor) (optional)
        @param c: Search depth for causal paths (optional) (default 3)
        @param p: Number of random permutations performed for significance analysis (optional) (default 1000)
        @param pagerank: Use Personalized PageRank to Diffuse (optional)
        @param all_paths: Use all paths instead of only causal paths (optional) (default False)
        @param singularity: if True, run using the Singularity container instead of the Docker container
        """

        if not edges or not sources or not targets or not output_file:
            raise ValueError("Required TieDIE arguments are missing")

        work_dir = "/spras"

        # Each volume is a tuple (src, dest) - data generated by Docker
        volumes = list()

        bind_path, edges_file = prepare_volume(edges, work_dir)
        volumes.append(bind_path)

        bind_path, sources_file = prepare_volume(sources, work_dir)
        volumes.append(bind_path)

        bind_path, targets_file = prepare_volume(targets, work_dir)
        volumes.append(bind_path)

        out_dir = Path(output_file).parent
        
        # TieDIE requires that the output directory exist
        out_dir.mkdir(parents=True, exist_ok=True)
        bind_path, mapped_out_dir = prepare_volume(str(out_dir), work_dir)
        volumes.append(bind_path) # Use posix path inside the container

        command = [
            "python",
            "/TieDIE/bin/tiedie",
            "-u",sources_file,
            "-d",targets_file,
            "-n",edges_file,
            "-s",str(s),
            "-c",str(c),
            "-p",str(p),
            "--pagerank", "True" if pagerank else "False",
            "--all_paths", "False" if all_paths else "False",
            "--output_folder", mapped_out_dir,
        ]

        print('Running TieDIE with arguments: {}'.format(' '.join(command)), flush=True)


        # TODO consider making this a string in the config file instead of a Boolean
        container_framework = "singularity" if singularity else "docker"
        out = run_container(container_framework, 
                            "reedcompbio/tiedie", 
                            command, 
                            volumes, 
                            work_dir)
        print(out)

        # Rename the primary output file to match the desired output filename
        output = Path(out_dir, "tiedie.sif")
        target = Path(output_file)
        output.rename(target)

    @staticmethod
    def parse_output(raw_pathway_file, standardized_pathway_file):
        """
        Convert a predicted pathway into the universal format
        @param raw_pathway_file: pathway file produced by an algorithm's run function
        @param standardized_pathway_file: the same pathway written in the universal format
        """
        print("Parsing TieDIE output")

        df = pd.read_csv(raw_pathway_file, sep="\t", names=[1, 2, 3])

        # get rid of the relationship column (since all relationships are the same : -a>)
        df_out = df.drop(columns=[2])
        df_out = add_rank_column(df_out)
        df_out.to_csv(standardized_pathway_file, sep="\t", index=False, header=False)
