"""
Description :

Function that performs the simplified persistence algorithm described in Sec.IV , resulting in a
list of birth - death pairs and the sector skeleton of the network .

Arguments :

edge_list : list of pairs of nodes comprising each edge
dP: list of absolute values of pressure difference for each edge
ascending : boolean value specifying ascending ( True ) or descending ( False ) filtration

Returns :

bd_pairs :list of birth - death pairs of pressure differences
bd_edges :list of birth - death pairs of edges
tau : list of persistence values for each birth - death pair
skeleton :set of edges that comprise sector skeleton

Note : Edges represented by positional index of edge in edge_list
"""
def persist_alg ( edge_list , dP, ascending = True ) :
    # Sort positional indices of edges in edge_list in ascending or descending order of absolute
    # value of pressure difference .
    NE = len( edge_list )
    edge_to_pressure= lambda edgei : dP[ edgei ]
    sorted_edges = sorted ( range ( NE ) , key = edge_to_pressure , reverse =( not ascending ) )
    
    # Map of sectors ( labeled by birth edges ) to sets of vertices seen so far in filtration .
    sectors_to_verts = dict()
    # Map of vertices seen so far in filtration to sectors .
    verts_to_sectors = dict()
    
    # Sector skeleton .
    skeleton = set()

    # Birth - death pairs of pressure differences .
    bd_pairs = []
    # Birth - death pairs of edges .
    bd_edges = []
    # Persistence values of each birth - death pair .
    tau = []

    # Iterate through each edge in sorted order .
    for i , ei in enumerate ( sorted_edges ) :
        ( vi , vj ) = edge_list [ ei ]
        # If both vertices have already appeared in filtration and are already in the same sector ,
        # then skip .
        if ( vi in verts_to_sectors and vj in verts_to_sectors
            and verts_to_sectors [ vi ] == verts_to_sectors [ vj ]) :
            continue
        # Add edge to sector skeleton .
        skeleton.add ( ei )
    
        # Case i ) : Edge creates new connected component .
        # Neither vertex has already appeared in filtration .
        # Create new sector and label sector with current edge as birth edge .
        if vi not in verts_to_sectors and vj not in verts_to_sectors :
            sectors_to_verts[ ei ] = { vi , vj }
            verts_to_sectors[ vi ] = ei
            verts_to_sectors[ vj ] = ei

        # Case ii ) : Edge merges two connected components .
        # Both vertices have already appeared in filtration and are in different sectors .
        # Record current edge as death edge , record persistence pair , and join sectors .
        elif ( vi in verts_to_sectors and vj in verts_to_sectors
            and verts_to_sectors [ vi ] != verts_to_sectors [ vj ]) :
        
            # Record current edge as death edge .
            death_edge = ei
            death_dP= dP[ ei ]
        
            si = verts_to_sectors[ vi ]
            sj = verts_to_sectors[ vj ]
        
            # Ensure that sector si is the more recent of the two sectors in the filtration .
            if ( ascending and dP[ si ] < dP[ sj ]) or ( not ascending and dP[ si ] > dP[ sj ]) :
                si , sj = sj , si
            
            # Record the label of this " younger " sector is the birth edge .
            # This is refered to as the " Elder Rule ."
            birth_edge = si
            birth_dP= dP[ si ]
        
            # Merge sector si , the younger sector , to sector sj , the older sector .
            merged_sector = sectors_to_verts.pop ( si )
            for vk in merged_sector :
                verts_to_sectors[ vk ] = sj
            
            sectors_to_verts[ sj ] = sectors_to_verts[ sj ] | merged_sector

            # Record persistence pair if persistence is greater that zero .
            if abs( death_dP- birth_dP) :
                bd_pairs.append(( birth_dP, death_dP) )
                bd_edges.append(( birth_edge , death_edge ) )
                tau.append( abs( death_dP- birth_dP) )

            # Case iii ) : Edge adds new vertex to the network .
            # Only one vertex has already appeared in filtration .
            # Add other vertex to same sector .
        else :#whats the problem?[]
	            if vi not in verts_to_sectors :
		            sj = verts_to_sectors[ vj ]
		            sectors_to_verts[ sj ].add ( vi )
		            verts_to_sectors[ vi ] = sj
	            else :
		            si = verts_to_sectors[ vi ]
		            sectors_to_verts[ si ].add ( vj )
		            verts_to_sectors[ vj ] = si
    return bd_pairs , bd_edges , tau , skeleton