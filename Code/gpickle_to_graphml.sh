#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-8:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p stats   # Partition to submit to
#SBATCH --mem=2000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors.err  # File to which STDERR will be written, %j inserts jobid
module load python/2.7.14-fasrc01
X=/n/regal/airoldi_lab/sussman/neurodata/
for i in `ls $X/gpickle`
do
    echo $X/$i | cut -d '.' -f 1
    python2 -c "import networkx as nx; g = nx.read_gpickle('`echo $X/gpickle/$i`'); nx.write_edgelist(g,'`echo $X/edgelist/$i | cut -d '.' -f 1`.edgelist',data=['weight'])"
    mv $X/gpickle/$i $X/gpickle_processed
    # nx.write_graphml(g,'`echo $X/graphml/$i | cut -d '.' -f 1`.graphml')"
done

# for i in `ls $1`
#	do
#	echo $1/graphml/$i | cut -d '.' -f 1
#	python2 -c "import networkx as nx; g = nx.read_gpickle('`echo $1/$i`'); nx.write_edgelist(g,'`echo $1/edgelist/$i| cut -d '.' -f 1`.edgelist',data=['weight'])"
    # nx.write_graphml(g,'`echo $1/graphml/$i | cut -d '.' -f 1`.graphml')"
# done

# import networkx as nx;
# fn = '/Volumes/Other/Data/neurodata_dtmri/MRN114sub-M87129789_ses-1_dwi_AAL.gpickle'
# g = nx.read_gpickle(fn);
# nx.write_graphml(g,'graphml/BNU1sub-0025864_ses-1_dwi_AAL.graphml')

# import networkx as nx;
# g = nx.read_gpickle("./SWU4_CPAC200_res-2x2x2_0025863_2.gpickle")
# nx.write_graphml(g,'graphml/SWU4_CPAC200_res-2x2x2_0025863_2.graphml')
# g.node

# write_edgelist(G, path, comments='#', delimiter=' ', data=True, encoding='utf-8')[source]
# Write graph as a list of edges.
# Parameters: 
# G : graph

# A NetworkX graph

# path : file or string

# File or filename to write. If a file is provided, it must be opened in ‘wb’ mode. Filenames ending in .gz or .bz2 will be compressed.

# comments : string, optional

# The character used to indicate the start of a comment

# delimiter : string, optional

# The string used to separate values. The default is whitespace.

# data : bool or list, optional

# If False write no edge data. If True write a string representation of the edge data dictionary.. If a list (or other iterable) is provided, write the keys specified in the list.

# encoding: string, optional

# Specify which encoding to use when writing file.
