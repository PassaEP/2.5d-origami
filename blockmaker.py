from collections import defaultdict
import math
import cadnano
from cadnano.document import Document
# from cadnano.proxies.cnenum import LatticeType  # GridType, StrandType
from cadnano.fileio.lattice import HoneycombDnaPart  # SquareDnaPart
from cadnano.part.nucleicacidpart import DEFAULT_RADIUS

from cadnano.part.refresholigoscmd import RefreshOligosCommand


# Honeycomb crossover locations
SCAF_LO_330_150 = [1, 11]  # \
SCAF_HI_330_150 = [2, 12]  # \
SCAF_LO_210_30 = [4, 15]   # /
SCAF_HI_210_30 = [5, 16]   # /
SCAF_LO_90_270 = [8, 18]   # |
SCAF_HI_90_270 = [9, 19]   # |

STAP_LO_330_150 = [6]  # \
STAP_HI_330_150 = [7]  # \
STAP_LO_210_30 = [20]  # /
STAP_HI_210_30 = [0]   # /
STAP_LO_90_270 = [13]  # |
STAP_HI_90_270 = [14]  # |

# Initial setup
radius = DEFAULT_RADIUS
doLattice = HoneycombDnaPart.legacyLatticeCoordToPositionXY
isEven = HoneycombDnaPart.isEvenParity

LO_PADDING = 21
HI_PADDING = 21
BUFFER_SPACE_Y = -6*radius  #down the coordinates become negative
BUFFER_SPACE_X = 2*(2*radius*math.sqrt(3)) #using the figma diagram
length_handles = 8

def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

def cleanCoords(coordList):
    cleanedList = []
    for coord in coordList:
        xyz = [truncate(coord[0], 3), truncate(coord[1], 3), truncate(coord[2], 3)]
        if xyz in cleanedList:
            continue

        cleanedList.append(xyz)

    return cleanedList

# Get a document
app = cadnano.app()

# takes a dictonary of design specs
# we need to generate the cadnano document
# ouputs



def cnFileMake(design):

    print(design['id'])

    doc = app.document = Document()

    # Create a new part
    part = doc.createNucleicAcidPart(use_undostack=False)

    # Manually track vh helix numbers
    vh_num_even = 0
    vh_num_odd = 1
    coord_to_vhnum = {}
    coord_to_vhnum_dummyhelices = {}
    vhnum_to_strandidxs = defaultdict(list)
    vhnum_to_strandidxs_dummy = defaultdict(list)
    vh_num_and_locations = {}

    # Create the VHs
    rows = design['rows']
    cols = design['cols']
    length = design['length']
    for row in range(rows):
        direction = 1 if row % 2 == 0 else -1
        for col in range(cols)[::direction]:
            # Determine xy coordinates
            x, y = doLattice(radius, row, col)
            x+=BUFFER_SPACE_X
            y+=BUFFER_SPACE_Y
            # Determine next vh_num
            if isEven(row, col):
                vh_num = vh_num_even
                vh_num_even += 2
            else:
                vh_num = vh_num_odd
                vh_num_odd += 2

            # Create the empty helix
            part.createVirtualHelix(x, y, 0., 21*(length//21)+42, id_num=vh_num)

            # Store vh_num for easy retrieval
            coord_to_vhnum[(row, col)] = vh_num

    initHiStrandIDs = defaultdict(list)
    initLoStrandIDs = defaultdict(list)
    vhnum_to_strandidxs = defaultdict(list)

    # Create initial strands in each row
    for row in range(rows):  # start a new loop to keep things clean
        direction = 1 if row % 2 == 0 else -1
        for col in range(cols)[::direction]:
            # Get the vh_num
            vh_num = coord_to_vhnum[(row, col)]

            # Get references to strandsets
            if isEven(row, col):
                scaf_strand_set, stap_strand_set = part.getStrandSets(vh_num)
            else:
                stap_strand_set, scaf_strand_set = part.getStrandSets(vh_num)

            # Calculate the low and high base idxs of the new scaf strand
            if direction == 1:
                # We are in an even row where vh_num is increasing left to right
                # and the even_vh <-> odd_ crossover angles will be 330° <-> 150°.
                # Also, the crossover will be on the HI idx of the crossover
                # since we are creating the "left" edge of the shape
                # We'll just default to the first location for now
                low_idx = SCAF_HI_330_150[0] + LO_PADDING
                # The high index is just the low index plus the length
                high_idx = SCAF_HI_330_150[0] + length - 1 + LO_PADDING
            else:
                # In odd rows, we are rastering back right to the left
                # so the xover locations for 210° <-> 30° angles will apply
                low_idx = SCAF_HI_210_30[0] + LO_PADDING
                high_idx = SCAF_HI_210_30[0] + length - 1 + LO_PADDING

            # Create the individual strands

            if type(high_idx) is float:
                high_idx = round(high_idx)

            if type(low_idx) is float:
                low_idx = round(low_idx)
                
            scaf_strand_set.createStrand(low_idx, high_idx, use_undostack=False)
            # print("Created strand: {}[{}...{}]".format(vh_num, low_idx, high_idx))

            vhnum_to_strandidxs[vh_num].append([low_idx, high_idx])
            stap_strand_set.createStrand(low_idx, high_idx, use_undostack=False)



    # Install edge crossovers between pairs of strands
    # currently only works for even col numbers

    for row in range(rows):  # new loop since we want strands in place first
        direction = 1 if row % 2 == 0 else -1
        for col in range(0, cols, 2)[::direction]:  # we will process pairs of strands
            # Get the vh_num
            if direction == 1:
                even_vh_num = coord_to_vhnum[(row, col)]
                odd_vh_num =  coord_to_vhnum[(row, col+1)]
            else:
                even_vh_num = coord_to_vhnum[(row, col+1)]
                odd_vh_num = coord_to_vhnum[(row, col)]

            # print('({},{}) {} {}'.format(row, col, even_vh_num, odd_vh_num))

            # To access strands we already created, we need to get the StrandSet
            # for each vh. We can discard the staple strandset for now
            even_scaf_strand_set, _ = part.getStrandSets(even_vh_num)
            _, odd_scaf_strand_set = part.getStrandSets(odd_vh_num)

            # Look up the the low and high indices of the strands we created
            even_lo_idx, even_hi_idx = vhnum_to_strandidxs[even_vh_num][0]
            odd_lo_idx, odd_hi_idx = vhnum_to_strandidxs[odd_vh_num][0]

            # This should match since we aimed to create equal-length strands above
            assert even_lo_idx == odd_lo_idx and even_hi_idx == odd_hi_idx

            # print('{}<->{} at {}={}, {}={}'.format(even_vh_num, odd_vh_num,
            #                                        even_lo_idx, odd_lo_idx,
            #                                        even_hi_idx, odd_hi_idx))

            # Get the reference to the strand we created on each helix
            even_strand = even_scaf_strand_set.getStrand(even_lo_idx)
            odd_strand = odd_scaf_strand_set.getStrand(odd_lo_idx)

            # Let's install the right (high) edge crossovers first
            e_idx3p = even_strand.idx3Prime()
            o_idx5p = odd_strand.idx5Prime()

            # We know based on the strand directions that
            # the even vh strand high idx should be its 3' end
            # the odd vh strand high idx should be its 5' end
            assert even_hi_idx == e_idx3p  # let's confirm anyway
            assert odd_hi_idx == o_idx5p

            # The 5' end of the even strand should align w/ the 3' of the odd strand
            assert e_idx3p == o_idx5p

            # Renaming these to make the directionality very clear
            from_strand, from_idx = even_strand, e_idx3p
            to_strand, to_idx = odd_strand, o_idx5p

            # Crossovers are ALWAYS created from the 3' end (triangle endpoint)
            # to the 5' end (rectangle endpoint).
            part.createXover(to_strand, to_idx,
                             from_strand, from_idx,
                             allow_reordering=True,
                             update_oligo=False,  # delayed for better performance
                             use_undostack=False)

            # Let's install the left (low) crossovers
            e_idx5p = even_strand.idx5Prime()
            o_idx3p = odd_strand.idx3Prime()

            # Renaming these to make the directionality very clear
            from_strand, from_idx = odd_strand, o_idx3p
            to_strand, to_idx = even_strand, e_idx5p

            part.createXover(to_strand, to_idx,
                             from_strand, from_idx,
                             allow_reordering=True,
                             update_oligo=False,  # delayed for better performance
                             use_undostack=False)

    # We need to call this after creating xovers with update_oligo=False
    RefreshOligosCommand(part).redo()



    # We could install the internal xovers by explicitly calculating the idxs
    # and handing the connections at direction reversal as a special edge case.
    # Let's install internal crossovers using a slightly different approach:
    # iterating over pairs of vhs, starting with the 2nd pair, ask cadnano
    # what possible crossover locations are available to connect back to the
    # # previous pair
    for vh_num in range(2, vh_num_even, 2):
        strand_idxs = vhnum_to_strandidxs[vh_num]
        print('INT-XO for', vh_num, strand_idxs)

        even_vh_num = vh_num
        odd_vh_num = vh_num-1

        even_scaf_strand_set, _ = part.getStrandSets(even_vh_num)
        _, odd_scaf_strand_set = part.getStrandSets(odd_vh_num)

        even_lo_idx, even_hi_idx = vhnum_to_strandidxs[even_vh_num][0]
        odd_lo_idx, odd_hi_idx = vhnum_to_strandidxs[odd_vh_num][0]

        approx_mid_idx = (even_lo_idx + even_hi_idx)//2
        print('even_id {}, prev_id: {} mid_idx: {}'.format(even_vh_num, odd_vh_num, approx_mid_idx))

        # potentialCrossoverMap returns:
        # Returns:
        #     dict: of :obj:`tuple` of form::
        #         neighbor_id_num: (fwd_hit_list, rev_hit_list)
        #     where each list has the form:
        #         [(id_num idx, [forward_neighbor_idxs], [reverse_neighbor_idxs]), ...]]

        per_neighbor_hits, pairs = part.potentialCrossoverMap(even_vh_num, approx_mid_idx)

        for neighbor_vh, hits in per_neighbor_hits.items():
            if neighbor_vh != odd_vh_num:
                continue

            fwd_hit_list, rev_hit_list = hits

            pxo_lo = fwd_hit_list[0]
            pxo_hi = fwd_hit_list[1]

            pxo_lo_from_idx = pxo_lo[0]
            pxo_lo_to_idx = pxo_lo[2][0]
            pxo_hi_from_idx = pxo_hi[2][0]
            pxo_hi_to_idx = pxo_hi[0]

            # Low idxs should match
            assert pxo_lo_from_idx == pxo_lo_to_idx
            # High idxs should match
            assert pxo_hi_from_idx == pxo_hi_to_idx

            # Low and High idxs should match
            assert pxo_hi_from_idx == pxo_lo_from_idx + 1
            assert pxo_hi_to_idx == pxo_lo_to_idx + 1

            # Split strands in preparation for crossovers
            # even_strand = even_scaf_strand_set.getStrand(pxo_lo_from_idx)
            # odd_strand = odd_scaf_strand_set.getStrand(pxo_lo_to_idx)
            # even_scaf_strand_set.splitStrand(even_strand, pxo_lo_from_idx)
            # odd_scaf_strand_set.splitStrand(odd_strand, pxo_hi_to_idx)

            # Low crossover
            from_strand = odd_scaf_strand_set.getStrand(pxo_lo_to_idx)
            from_idx = pxo_lo_from_idx
            to_strand = even_scaf_strand_set.getStrand(pxo_lo_from_idx)
            to_idx = pxo_lo_to_idx
            part.createXover(to_strand, to_idx,
                             from_strand, from_idx,
                             allow_reordering=False,
                             update_oligo=False,  # delayed for better performance
                             use_undostack=False)

            # High crossover
            from_strand = even_scaf_strand_set.getStrand(pxo_hi_from_idx)
            from_idx = pxo_hi_from_idx
            to_strand = odd_scaf_strand_set.getStrand(pxo_hi_to_idx)
            to_idx = pxo_hi_to_idx
            part.createXover(to_strand, to_idx,
                             from_strand, from_idx,
                             allow_reordering=False,
                             update_oligo=False,  # delayed for better performance
                             use_undostack=False)

        # We need to call this after creating xovers with update_oligo=False
        RefreshOligosCommand(part).redo()

    # Autostaple
    for vh_num in range(0, vh_num_even, 2):
        strand_idxs = vhnum_to_strandidxs[vh_num]
        stap_lo_idx, stap_hi_idx = strand_idxs[0]
        # print('Installing staple XOs for', vh_num, strand_idxs)

        # even_lo_idx, even_hi_idx = vhnum_to_strandidxs[even_vh_num][0]
        # odd_lo_idx, odd_hi_idx = vhnum_to_strandidxs[odd_vh_num][0]
        # get a list of n
        per_neighbor_hits, pairs = part.potentialCrossoverMap(vh_num)

        # Since we are always installing crossovers from an even vh,
        # the staple is always on the reverse strand (5'-3' right to left)

        for neighbor_vh, hits in per_neighbor_hits.items():
            _, rev_hit_list = hits

            # hit[1] is the forward strand of the neighbor
            all_staple_xover_idxs = [hit[1][0] for hit in rev_hit_list]

            staple_xover_idxs = list(
                filter(
                    lambda x: x > stap_lo_idx+6 and x < stap_hi_idx-6,
                    all_staple_xover_idxs)
            )
            # print("  {} is a neighbor".format(neighbor_vh))
            # print("  filtered staple hits", staple_xover_idxs)

            stap_idx = 0
            while stap_idx < len(staple_xover_idxs)-1:
                lo_idx = staple_xover_idxs[stap_idx]
                hi_idx = staple_xover_idxs[stap_idx+1]
                if hi_idx - lo_idx != 1:
                    stap_idx += 1
                    continue

                # print(lo_idx, hi_idx, end=' ')

                # even local ss
                _, vh_stap_ss = part.getStrandSets(vh_num)
                # odd neighbor ss
                neighbor_stap_ss, _ = part.getStrandSets(neighbor_vh)

                # Low crossover
                from_strand = vh_stap_ss.getStrand(lo_idx)
                from_idx = lo_idx
                to_strand = neighbor_stap_ss.getStrand(lo_idx)
                to_idx = lo_idx
                part.createXover(to_strand, to_idx,
                                 from_strand, from_idx,
                                 allow_reordering=False,
                                 update_oligo=False,  # delayed for better performance
                                 use_undostack=False)

                # # High crossover
                from_strand = neighbor_stap_ss.getStrand(hi_idx)
                from_idx = hi_idx
                to_strand = vh_stap_ss.getStrand(hi_idx)
                to_idx = hi_idx
                part.createXover(to_strand, to_idx,
                                 from_strand, from_idx,
                                 allow_reordering=False,
                                 update_oligo=False,  # delayed for better performance
                                 use_undostack=False)

                stap_idx += 2

            RefreshOligosCommand(part).redo()


    #DUMMY STUFF STARTS HERE

    #Record virtual helix number of the dummies that we want

    dummyVH = []

    hiDummyXO = []

    loDummyXO = []

    # Add dummy helices : first row + 2r, last row -2r

    for row in range(rows):
        direction = 1 if row % 2 == 0 else -1
        if row == 0:  # upper dummy helices
            for col in range(cols)[::direction]:
                if col % 2 == 0:
                    x, y = doLattice(radius, row, col)
                    x += BUFFER_SPACE_X
                    y += BUFFER_SPACE_Y + 2*radius
                # Determine next vh_num
                    vh_num = vh_num_odd
                    vh_num_odd += 2

                    # Create the empty helix
                    part.createVirtualHelix(x, y, 0., 21*(length//21)+42, id_num=vh_num)

                # Store vh_num for easy retrieval
                    coord_to_vhnum_dummyhelices[(row, col)] = vh_num
                    dummyVH.append(vh_num)

        if row == rows-1:
            if direction == 1:
                for col in range(cols)[::direction]:
                    if col % 2 == 1:
                        x, y = doLattice(radius, row, col)
                        x += BUFFER_SPACE_X
                        y += BUFFER_SPACE_Y - 2*radius
                    # Determine next vh_num
                        vh_num = vh_num_even
                        vh_num_even += 2
                        part.createVirtualHelix(x, y, 0., 21*(length//21)+42, id_num=vh_num)
                        coord_to_vhnum_dummyhelices[(row, col)] = vh_num
                        print(coord_to_vhnum_dummyhelices,row,col)
                        dummyVH.append(vh_num)

            else:
                for col in range(cols)[::direction]:
                    if col % 2 == 0:
                        x, y = doLattice(radius, row, col)
                        x += BUFFER_SPACE_X
                        y += BUFFER_SPACE_Y - 2*radius
                        vh_num = vh_num_even
                        vh_num_even += 2
                        part.createVirtualHelix(x, y, 0., 21*(length//21)+42, id_num=vh_num)
                        coord_to_vhnum_dummyhelices[(row, col)] = vh_num
                        dummyVH.append(vh_num)


                # Create the empty helix

    # For the dummy helices
    staple_xover_idxs = []
    for vh_num in range(cols*rows, vh_num_even, 1):
        per_neighbor_hits, pairs = part.potentialCrossoverMap(vh_num)
        for neighbor_vh, hits in per_neighbor_hits.items():

            # Need to match [even|odd] dummy helices with [rev|fwd] staple ss
            if vh_num % 2 == 0:
                _, staple_hit_list = hits
                neighbor_stap_strand = 1
            else:
                staple_hit_list, _ = hits
                neighbor_stap_strand = 2

            # hit[1] is the forward strand of the neighbor
            all_staple_xover_idxs = [hit[neighbor_stap_strand][0] for hit in staple_hit_list]

            staple_xover_idxs = list(
                filter(
                    lambda x: x > stap_lo_idx+6 and x < stap_hi_idx-6,
                    all_staple_xover_idxs)
            )

            # We again need to filter out "lone" staple xovers -SD
            paired_staple_xo_idxs = []
            i = 0
            while i < len(staple_xover_idxs)-1:
                lo_idx = staple_xover_idxs[i]
                hi_idx = staple_xover_idxs[i+1]
                if hi_idx - lo_idx != 1:
                    i += 1
                    continue
                paired_staple_xo_idxs.extend([lo_idx, hi_idx])
                i += 2

            print("  {} has neighbor {}".format(vh_num, neighbor_vh))
            print("  filtered staple hits", paired_staple_xo_idxs)
            vh_num_and_locations[vh_num] = paired_staple_xo_idxs

    for row in [0, rows-1]:
        direction = 1 if row % 2 == 0 else -1
        for col in range(cols)[::direction]:
            if row == 0:
                if col % 2 == 0:
                    # Get the vh_num
                    vh_num = coord_to_vhnum_dummyhelices[(row, col)]
                    # Get references to strandsets
                    if isEven(row, col):
                        # _, stap_strand_set = part.getStrandSets(vh_num)
                        stap_strand_set, _ = part.getStrandSets(vh_num)
                    else:
                        # stap_strand_set, _ = part.getStrandSets(vh_num)
                        _, stap_strand_set = part.getStrandSets(vh_num)
                    # Calculate the low and high base idxs of the new scaf strand
                    if direction == 1:
                        # We are in an even row where vh_num is increasing left to right
                        # and the even_vh <-> odd_ crossover angles will be 330° <-> 150°.
                        # Also, the crossover will be on the HI idx of the crossover
                        # since we are creating the "left" edge of the shape
                        # We'll just default to the first location for now

                        # The high index is just the low index plus the length
                        for i in range(1, len(vh_num_and_locations[vh_num])-1, 2):
                            low_idx = vh_num_and_locations[vh_num][i]
                            high_idx = low_idx+length_handles
                            vhnum_to_strandidxs_dummy[vh_num].append((low_idx, high_idx))
                            print("creating handle {}[{}]-[{}]".format(vh_num, low_idx, high_idx))
                            stap_strand_set.createStrand(low_idx, high_idx, use_undostack=False)
                print("test22",stap_strand_set.getStrand(35))
            elif row == rows-1:
                if rows % 2 != 0:
                    if col % 2 != 0:
                        # Get the vh_num
                        vh_num = coord_to_vhnum_dummyhelices[(row, col)]
                        # Get references to strandsets
                        if isEven(row, col):
                            # _, stap_strand_set = part.getStrandSets(vh_num)
                            stap_strand_set, _ = part.getStrandSets(vh_num)
                        else:
                            # stap_strand_set, _ = part.getStrandSets(vh_num)
                            _, stap_strand_set = part.getStrandSets(vh_num)
                        # Calculate the low and high base idxs of the new scaf strand
                        if direction == 1:
                            for i in range(0, len(vh_num_and_locations[vh_num])-1, 2):
                                high_idx = vh_num_and_locations[vh_num][i]
                                low_idx = high_idx-length_handles
                                vhnum_to_strandidxs_dummy[vh_num].append((low_idx, high_idx))
                                print("creating handle {}[{}]-[{}]".format(vh_num, low_idx, high_idx))
                                stap_strand_set.createStrand(low_idx, high_idx, use_undostack=False)
                else:
                    if col % 2 == 0:
                        print(row,col)
                    # Get the vh_num
                        vh_num = coord_to_vhnum_dummyhelices[(row, col)]
                        # Get references to strandsets
                        if isEven(row, col):
                            # _, stap_strand_set = part.getStrandSets(vh_num)
                            stap_strand_set, _ = part.getStrandSets(vh_num)
                        else:
                            # stap_strand_set, _ = part.getStrandSets(vh_num)
                            _, stap_strand_set = part.getStrandSets(vh_num)
                        # Calculate the low and high base idxs of the new scaf strand
                        if direction == -1:
                            for i in range(0,len(vh_num_and_locations[vh_num])-1, 2):
                                high_idx = vh_num_and_locations[vh_num][i]
                                low_idx = high_idx-length_handles
                                vhnum_to_strandidxs_dummy[vh_num].append((low_idx, high_idx))
                                print("creating handle {}[{}]-[{}]".format(vh_num, low_idx, high_idx))
                                stap_strand_set.createStrand(low_idx, high_idx, use_undostack=False)
                                #installing crossovers
    staple_xover_idxs = []
    idx_vhnum = {}
    it = 0

    xoBox = defaultdict(list)

    for vh_num in range(cols*rows, vh_num_even, 1):
        print("vh_num",vh_num)
        per_neighbor_hits, pairs = part.potentialCrossoverMap(vh_num)
        for neighbor_vh, hits in per_neighbor_hits.items():

            # Need to match [even|odd] dummy helices with [rev|fwd] staple ss
            if vh_num % 2 == 0:
                _, staple_hit_list = hits
                neighbor_stap_strand = 1
            else:
                staple_hit_list, _ = hits
                neighbor_stap_strand = 2

            # hit[1] is the forward strand of the neighbor
            all_staple_xover_idxs = [hit[neighbor_stap_strand][0] for hit in staple_hit_list]

            staple_xover_idxs = list(
                filter(
                    lambda x: x > stap_lo_idx+6 and x < stap_hi_idx-6,
                    all_staple_xover_idxs)
            )

            # We again need to filter out "lone" staple xovers -SD
            paired_staple_xo_idxs = []
            i = 0
            while i < len(staple_xover_idxs)-1:
                lo_idx = staple_xover_idxs[i]
                hi_idx = staple_xover_idxs[i+1]
                if hi_idx - lo_idx != 1:
                    i += 1
                    continue
                paired_staple_xo_idxs.extend([lo_idx, hi_idx])
                i += 2

            print("  {} has neighbor {}".format(vh_num, neighbor_vh))
            print("  filtered staple hits", paired_staple_xo_idxs)
            vh_num_and_locations[vh_num] = paired_staple_xo_idxs

            # stap_strand_set, _ = part.getStrandSets(vh_num)
            if vh_num % 2 !=0:
                # _, stap_strand_set = part.getStrandSets(vh_num)
                stap_strand_set, _ = part.getStrandSets(vh_num)
                print("uses iseven")
            else:
                # stap_strand_set, _ = part.getStrandSets(vh_num)
                _, stap_strand_set = part.getStrandSets(vh_num)
                print("else")
            # print(part.getStrandSets(5))
            # print(vh_num,"part",row, col,stap_strand_set.getStrand())
            # # even local ss
            if vh_num % 2 == 0:
            # odd neighbor ss
                neighbor_stap_ss, _ = part.getStrandSets(neighbor_vh)

                for i in range(0,len(vh_num_and_locations[vh_num])-1, 2):
                    high_idx = vh_num_and_locations[vh_num][i]
                    from_strand = stap_strand_set.getStrand(high_idx)
                    from_idx = high_idx
                    to_strand = neighbor_stap_ss.getStrand(high_idx)
                    to_idx = high_idx


                    print("creating xover from {} to {} at {} = {}".format(from_strand,to_strand,from_idx, to_idx))

                    part.createXover(to_strand, to_idx,
                                     from_strand, from_idx,
                                     allow_reordering=False,
                                     update_oligo=False,  # delayed for better performance
                                     use_undostack=False)
                    print("5' prime end", from_idx,vh_num) #5' end
                    idx_vhnum[vh_num] = from_idx
                    xoBox[vh_num].append(from_idx)

            else:
                print(vh_num,"vhnum")
                _, neighbor_stap_ss = part.getStrandSets(neighbor_vh)

                for i in range(1,len(vh_num_and_locations[vh_num])-1, 2):
                    low_idx = vh_num_and_locations[vh_num][i]
                    to_strand = neighbor_stap_ss.getStrand(low_idx)
                    from_idx = low_idx
                    from_strand = stap_strand_set.getStrand(low_idx)
                    to_idx = low_idx
                    print("vhnum",vh_num)
                    print("creating xover from {} to {} at {} = {}".format(from_strand,to_strand,from_idx, to_idx))
                    part.createXover(to_strand, to_idx,
                                     from_strand, from_idx,
                                     allow_reordering=False,
                                     update_oligo=False,  # delayed for better performance
                                     use_undostack=False)
                    print("5' prime end", from_idx,vh_num)
                    idx_vhnum[vh_num] = from_idx
                    xoBox[vh_num].append(from_idx)



        RefreshOligosCommand(part).redo()

        handleEndCoords = []
        for dummyNum in dummyVH:
            # traverse through the list
            ids = xoBox[dummyNum]
            if dummyNum % 2 == 0: # if even, get reverse
                for id in ids:
                    _, _, revStapleStrands = part.getCoordinates(dummyNum)
                    handleEndCoords.append(revStapleStrands[id])
            else:
                for id in ids:
                    _, fwdStapleStrands, _ = part.getCoordinates(dummyNum)
                    handleEndCoords.append(fwdStapleStrands[id])


        # cleanedHandleCoords = cleanCoords(handleEndCoords)

        print('different counts of crossovers')
        print(xoBox)
        print(dummyVH)

        # print('minimum and max xy coords')
        #
        minX, minY, maxX, maxY = part.boundDimensions()
        #
        # print('handle coords')
        #
        # print(cleanedHandleCoords)

        zList = []

        for vh_num in range(0, cols*rows):
            _, set2, set3 = part.getCoordinates(vh_num)
            for coordinate in set2:
                zList.append(float(coordinate[2]))
            for coordinate in set3:

                zList.append(float(coordinate[2]))



        minZ = min(zList)

        maxZ = max(zList)


        return minX, minY, minZ, maxX, maxY, maxZ, handleEndCoords


    doc.writeToFile('scripted_design_{}.json'.format(design['id']), legacy=True)



if __name__ == '__main__':
    sampleDesign = {'id': 1, 'rows': 2, 'cols': 3, 'length': 84}
    cnFileMake(sampleDesign)
