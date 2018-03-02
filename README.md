# meta3Cbox 

A set of scripts that streamlines the processing and binning of metagenomic 3C datasets.
    
Written by L. Baudry, T. Foutel-Rodier and M. Marbouty (with contributions from A. Cournac and V. Scolari)

*Spatial Regulation of Genomes* (Institut Pasteur, Paris)
    
Version 0.1a [March 2018]
    
## Usage


    ./meta3c.sh {align|partition|annotation|binning|pipeline} [parameters]

    
A meta3Cbox command takes the form './meta3c.sh action --param1 arg1 --param2 arg2' etc.
    
There are four actions or steps in the meta3c pipeline. They must be run in this order:
    
* align: map paired-end reads on a preliminary assembly, then generate a network from
 detected contacts between DNA chunks.
* partition: perform the Louvain community detection algorithm many times to isolate
     chunks that consistently cluster together for binning purposes.
* annotation: run standard annotation software on the assembly (namely gene prediction
      and database comparison) to match with the bins
* binning: match annotations to bins, extract bin genomes and subnetworks, build bin-local
   and global contact maps
    
After the last step is completed, you should have at your disposal a set of bins, their relative
enrichments in various gene categories as well as the contact map of each bin.
    
In addition, there are a number of optional, miscellaneous actions:
    
* pipeline: check the environment is right, then run all four of the above sequentially.
    This can take a while.
* deploy: set up the environment on Ubuntu 14.04 and higher
* version: display current version number
* help: display this (hopefully useful) help message
    
   Please refer to the meta3Cbox manual for detailed explanations on the parameters.
    
## Contact

lyam.baudry@pasteur.fr or romain.koszul@pasteur.fr
    

