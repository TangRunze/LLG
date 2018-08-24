
for i in `ls $1`
	do
	echo $1/graphml/$i | cut -d '.' -f 1
	python2 -c "import networkx as nx; g = nx.read_gpickle('`echo $1/$i`'); nx.write_edgelist(g,'`echo $1/edgelist/$i | cut -d '.' -f 1`.edgelist',data=['weight'])"
    # nx.write_graphml(g,'`echo $1/graphml/$i | cut -d '.' -f 1`.graphml')"
done

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