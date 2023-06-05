import pandas as pd
import warnings
from src.prm import PRM
from pathlib import Path
from src.util import prepare_volume, run_container

__all__ = ['RandomWalk']

class RandomWalk(PRM):
    required_inputs = ['edges', 'sources', 'targets']

    @staticmethod
    def generate_inputs(data, filename_map):
        """
        Access fields from the dataset and write the required input files
        @param data: dataset
        @param filename_map: a dict mapping file types in the required_inputs to the filename for that type
        @return:
        """
        # ensures the required input are within the filename_map 
        for input_type in RandomWalk.required_inputs:
            if input_type not in filename_map:
                raise ValueError(f"{input_type} filename is missing")

        # will take the sources and write them to files, and repeats with targets 
        for node_type in ['sources', 'targets']:                
            nodes = data.request_node_columns([node_type])  
            if data.contains_node_columns('prize'):
                #NODEID is always included in the node table
                node_df = data.request_node_columns(['prize'])
                nodes = pd.merge(nodes, node_df, on='NODEID')
                # creates with the node type without headers 
                nodes.to_csv(filename_map[node_type], index=False, sep= " ", columns=['NODEID', 'prize'])
            else:
                #If there aren't prizes but are sources and targets, make prizes based on them
                nodes = data.request_node_columns([node_type])
                # make all nodes have a prize of 1
                nodes['prize'] = 1.0
                # creates with the node type without headers 
                nodes.to_csv(filename_map[node_type], index=False, sep= " ", columns=['NODEID', 'prize'])

        # create the network of edges 
        edges = data.get_interactome()
        
        # creates the edges files that contains the head and tail nodes and the weights after them
        edges.to_csv(filename_map['edges'], sep=" ", index=False, columns=["Interactor1","Interactor2","Weight"])
        
        
    # Skips parameter validation step
    @staticmethod
    def run(edges=None, sources=None, targets = None, output_file = None, df = 0.85, singularity=False):
        """
        Run LocalNeighborhood with Docker
        @param nodetypes:  input node types with sources and targets (required)
        @param network:  input network file (required)
        @param output_file: path to the output pathway file (required)
        @param k: path length (optional)
        @param singularity: if True, run using the Singularity container instead of the Docker container
        """
       
        if not edges or not sources or not targets or not output_file:
            raise ValueError('Required RandomWalk arguments are missing')

        work_dir = '/spras'

        # Each volume is a tuple (src, dest) - data generated by Docker
        volumes = list()
  
        bind_path, edges_file = prepare_volume(edges, work_dir)
        volumes.append(bind_path)

        bind_path, sources_file = prepare_volume(sources, work_dir)
        volumes.append(bind_path)
        
        bind_path, targets_file = prepare_volume(targets, work_dir)
        volumes.append(bind_path)

        print("REACHED HERE")

        # PathLinker does not provide an argument to set the output directory
        # Use its --output argument to set the output file prefix to specify an absolute path and prefix
        out_dir = Path(output_file).parent
        # PathLinker requires that the output directory exist
        out_dir.mkdir(parents=True, exist_ok=True)
        bind_path, mapped_out_dir = prepare_volume(str(out_dir), work_dir)
        volumes.append(bind_path)
        mapped_out_prefix= mapped_out_dir + '/out'  # Use posix path inside the container
        

        # need 2 output files, one for nodes and one for edges (not implemented yet)
        # if inputs are changed to output_nodes and output_edges, then the error will be "unexpected argument : output_file"
        # if inputs are changed to only output_file, then the final output file will contain edges.
        command = ['python',
                   '/RandomWalk/random_walk.py',
                   '--edges_file', edges_file,
                   '--sources_file', sources_file,
                   '--targets_file', targets_file,
                   '--output_prefix', mapped_out_prefix]

        print('Running RandomWalk with arguments: {}'.format(' '.join(command)), flush=True)

        print("HERE")
        # TODO consider making this a string in the config file instead of a Boolean
        container_framework = 'singularity' if singularity else 'docker'
        out = run_container(container_framework,
                            'erikliu24/random-walk',
                            command,
                            volumes,
                            work_dir)
        print(out)

        # Rename the primary output file to match the desired output filename
        # Currently PathLinker only writes one output file so we do not need to delete others
        # We may not know the value of k that was used
        output = Path(out_dir, 'out')
        output.rename(output_file)
        
    @staticmethod
    def parse_output(raw_pathway_file, standardized_pathway_file):
        """
        Convert a predicted pathway into the universal format
        @param raw_pathway_file: pathway file produced by an algorithm's run function
        @param standardized_pathway_file: the same pathway written in the universal format
        """
        print('Parsing LocalNeighborhood output')
        df = pd.read_csv(raw_pathway_file, sep='|', header=None)
        df.to_csv(standardized_pathway_file, header=False, index=False, sep='|')
