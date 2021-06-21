# (C) 2011, Alexander Korsunsky
# (C) 2011-2015 Elke Schaper

"""
    :synopsis: Parsing repeat detection algorithm output

    .. moduleauthor:: Elke Schaper <elke.schaper@isb-sib.ch>
"""

import logging
import re
import csv
from collections import OrderedDict

LOG = logging.getLogger(__name__)


class RepeatRegion:

    def __init__(self, protein_id="", begin=None, msa=None):
        self.protein_id = protein_id
        self.begin = begin  # 1-based
        if msa is None:
            msa = []
        self.msa = msa


def tred_get_repeats(infile):
    r"""Read repeats from a TRED standard output (stdout) file stream successively.

    Read repeats from a TRED standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of TRED output file::

         Start: start End: \d+ Length: \d+

        ( \d repeat_unit \d
        ( alignment_indicator )?)*

    Args:
        infile (file stream): File stream of output1 from
            tred1 myDNA.faa intermediate_output
            tred2 myDNA.faa intermediate_output output1 output2

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout TRED output syntax.
    """

    # pattern for start index
    pat_start = re.compile(r"Start: (\d+)")
    pat_repeat_unit = re.compile(r"\d+\s+([ACGT-]+)")
    pat_alignment_indicator = re.compile(r"\s+([ACGT-]+)")

    # Our possible parser states:
    #
    # 1: searching for start
    # 2: searching for first repeat unit
    # 3: searching for all further repeat units

    identifier = ''
    state = 1
    repeat_units = []
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])
        if 1 == state:
            region = RepeatRegion(protein_id=identifier, msa=[])
            match = pat_start.match(line)
            if match:
                LOG.debug(" * (1->2) Found start")
                LOG.debug("Start: %s", match.group(1))
                region.begin = int(match.group(1))
                state = 2

        if 2 == state:
            match = pat_repeat_unit.match(line)
            if match:
                LOG.debug(" * (2->3) Found first repeat_unit")
                repeat_units.append((match.group(1), True))
                state = 3
                continue

        if 3 == state:
            match = pat_repeat_unit.match(line)
            if match:
                LOG.debug(" * (3->3) Found another repeat_unit")
                repeat_units.append((match.group(1), True))
            else:
                match = pat_alignment_indicator.match(line)
                if match:
                    LOG.debug(" * (3->3) Found an alignment_indicator unit")
                    repeat_units.append((match.group(1), False))
                else:
                    LOG.debug(" * (3->1) Found end of repeat (yielding)")
                    state = 1
                    region.msa = tred_msa_from_pairwise(repeat_units)
                    yield region


def tred_msa_from_pairwise(repeat_units):
    """ Construct a MSA from pairwise alignments.

    Construct a MSA from pairwise alignments. At the moment, gaps following the repeat are
    not added. However, these gaps are added automatically when a ``Repeat`` instance is
    created.

    Args:
        repeat_units (list of str): Read in from TRED output files

    Returns:
         (list of str)

    .. todo:: Is the Args format mentioned correctly?
    """

    pat_gap = re.compile(r"(-*)")
    result = []
    index = 0
    for iR in range(len(repeat_units)):
        ru = repeat_units[iR]

        # The next repeat unit
        if ru[1]:
            result.append('-' * index + ru[0])
            # How many gaps in the beginning of this repeat unit?
            index += len(pat_gap.match(ru[0]).group())

        # The next alignment indicator
        else:
            for iL in range(len(ru[0])):
                if ru[0][iL] == '-':
                    # enter a gap between
                    # the index + iL and index + iL + 1 character
                    # in each repeat unit in result so far:
                    result = [iRU[:index + iL] + '-' + iRU[index + iL:]
                              for iRU in result]

    return result


def treks_get_repeats(infile):
    """ Read repeats from a T-REKS standard output (stdout) file stream successively.

    Read repeats from a T-REKS standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of T-REKS output file::

        protein ::=
            ">" identifier
            repeat*
        #
        repeat ::=
            repeat_header
            sequence*
            "*"+
        #
        repeat_header ::= "Length:" integer "residues - nb:" integer  "from"  integer "to" integer "- Psim:"float "region Length:"integer


    Args:
        infile (file stream): File stream to the file of the standard output of T-Reks

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout T-REKS output syntax.
    """

    # pattern for protein identifier
    pat_identifier = re.compile(r">(\S+)")

    # pattern for repeat properties
    pat_repeat_header = re.compile(
        r"Length:\s*\d+\s+"
        r"residues - nb:\s*(\d+)\s+"  # 1. number of repeats
        r"from\s+(\d+)\s+to\s+(\d+)\s*"  # 2-3. start and end of repeat region (1-based)
        r"-\s*Psim:\s*(-?[\d\.]+)\s+"  # 4. probability (accept negative for (malformed) log10 values)
        r"region Length:\s*(\d+)")  # 5. length

    pat_repeat_end = re.compile(r"\*+")

    # pattern for repeat sequence
    # FIXME stuff that occurs here is not just A-Z but a set of possible amino
    # acid symbols
    pat_sequence = re.compile(r"([A-Z-]+)")

    # Our possible parser states:
    #
    # 1: state between repeats
    #   entry: reset repeat
    #   expect repeat_header(goto 2) OR identifier(store identifier, goto 1)
    # 2: state for multiple sequence alignment line
    # expect sequence(goto 3, append sequence) OR repeat_end(return repeat,
    # goto 1)

    state = 1
    identifier = ""

    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])

        if 1 == state:  # Parsing header
            region = RepeatRegion(protein_id=identifier, msa=[])
            LOG.debug("msa: %s", "\n".join(region.msa))
            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" * (1->2) Found properties")
                region.begin = int(match.group(2))  # 1-based in T-Reks output
                end = int(match.group(3))
                state = 2

            match = pat_identifier.match(line)
            if match:
                LOG.debug(" * (1->1) Found identifier")
                identifier = match.group(1)
                state = 1

        elif 2 == state:  # Parsing MSA
            match = pat_sequence.match(line)
            if match:
                LOG.debug(" * (2->2) Found MSA line (appending)")
                region.msa.append(match.group(1))
                state = 2
            match = pat_repeat_end.match(line)
            if match:
                LOG.debug(" * (2->1) Found end of repeat (yielding)")
                state = 1

                # Perform some checks before yielding
                msalen = sum(1 for row in region.msa for c in row if c != '-')
                if region.begin + msalen != end + 1:
                    logging.warn("Inconsistency in T-Reks file: MSA length doesnt match start and end")

                yield region

        else:
            raise AssertionError("Huh? Unknown parser state " +
                                 str(state))


def xstream_get_repeats(infile):
    """ Read repeats from a XSTREAM output xls chart

    Read repeats from a XSTREAM output xls chart

    Postcondition: infile points to EOF.

    Args:
        infile (file stream): File stream to read XSTREAM output from a xls chart.

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.
    """

    # The infile is luckily enough in csv format:
    reader = csv.reader(infile, dialect='excel-tab')

    # read first row with the fieldnames
    header = next(reader)

    LOG.debug("Header has %d fields.", len(header))

    if len(header) != 8 and len(header) != 10:
        raise ValueError(
            "XStream output file seems to be malformed. "
            "Make sure to open this function with the generated xls chart.")

    # when XSTREAM is fed files with multiple sequences, first field is "identifier",
    # otherwise it is "seq length"
    if header[0] == "identifier":
        field_offset = 2
    else:
        field_offset = 0

    for row in reader:
        region = RepeatRegion()
        if field_offset:
            region.protein_id = row[0]

        region.begin = int(row[0 + field_offset])
        region.msa = row[4 + field_offset].split()
        LOG.debug("Found repeat, yielding.")
        if len(region.msa) >= 2:
            yield region


def trust_fill_repeats(msa, begin, sequence, maximal_gap_length=20):
    ''' return a trust msa that has no longer indels than maximal_gap_length,
    that contains the indel characters even when not part of the trust output file.
    Background trust returns tandem repeats, but also distant repeats. '''
    gapless_msa = [repeat_unit.replace('-', '').upper() for repeat_unit in msa]
    sequence = sequence.upper()

    # Find the start and end positions of the predicted repeat units
    position = [(begin - 1, begin + len(gapless_msa[0]) - 2)]
    for repeat_unit in gapless_msa[1:]:
        find_index = sequence[position[-1][1] + 1:].find(repeat_unit)
        # Repeat unit could not be found in sequence -> Discard the repeat
        if find_index == -1:
            return None, None
        repeat_unit_begin = find_index + position[-1][1] + 1
        position.append(
            (repeat_unit_begin, repeat_unit_begin + len(repeat_unit) - 1))

    # Derive the start and end positions of the gaps
    gap_position = [(i[1] + 1, j[0] - 1)
                    for i, j in zip(position[:-1], position[1:])]
    gaps = [i[1] - i[0] + 1 for i in gap_position]

    # Filter out repeat units that are further apart than maximal_gap_length
    gap_valid = ''.join(
        ['1' if i_gap <= maximal_gap_length else '0' for i_gap in gaps])

    count_valid_pairs = [len(m.group())
                         for m in re.finditer(re.compile(r'1+'), gap_valid)]
    # All repeat units are further apart than maximal_gap_length? -> Discard
    # the repeat
    if len(count_valid_pairs) == 0:
        return None, None

    # Choose the sequence of pairs closer then maximal_gap_length that is
    # longest
    valid_index = gap_valid.find('1' * max(count_valid_pairs))

    # Shorten the predicted msa accordingly.
    msa = msa[valid_index:valid_index + max(count_valid_pairs) + 1]
    gaps = gaps[valid_index:valid_index + max(count_valid_pairs) + 1]
    gap_position = gap_position[
        valid_index:valid_index +
        max(count_valid_pairs) +
        1]
    position = position[valid_index:valid_index + max(count_valid_pairs) + 1]

    # Add missing sequence to the repeat units
    gap_count_before = 0
    for i, i_gap in enumerate(gaps):
        gap_count_after = gap_count_before + i_gap
        msa[i] += '-' * (gap_count_before) + sequence[gap_position[i][0]:gap_position[i][1] + 1] + '-' * (sum(gaps) - gap_count_after)
        gap_count_before = gap_count_after
    msa[-1] += '-' * sum(gaps)

    return msa, position[0][0] + 1


def trust_get_repeats(infile):
    """ Read repeats from a TRUST standard output (stdout) file stream successively.

    Read repeats from a TRUST standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of TRUST standard output::

        protein ::=
            ">" identifier
            (repeat_types)*
            "//"
        #
        repeat_types ::=
            "REPEAT_TYPE" integer
            "REPEAT_LENGTH" integer
            (repeat_info)*
            (">Repeat " integer
            sequence)*
        #
        repeat_info ::=
            integer integer [integer] [integer]

    Args:
        infile (file stream): File stream from TRUST standard output.

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout TRUST output syntax.
    """

    pat_identifier = re.compile(r">(\S+)")
    pat_repeat_type = re.compile(r"REPEAT_TYPE \d+")
    pat_repeat_length = re.compile(r"REPEAT_LENGTH \d+")
    pat_repeat_info = re.compile(r"(\d+) (\d+).*$")
    pat_repeat_header = re.compile(r">Repeat \d+")
    # FIXME find proper character set here
    pat_repeat_sequence = re.compile(r"([A-Za-z-]+)")
    pat_protein_end = re.compile(r"//")

    # Our possible parser states:
    #
    # 0: initial state
    #   expect: identifier(store identifier, goto 1)
    # 1: before the beginning of a repeat
    #   expect: "REPEAT_TYPE"(goto 2) or "//"(goto 0)
    # 2: state after having found "REPEAT_TYPE"
    #   entry: reset state of return value (repeat)
    #   expect: "REPEAT_LENGTH"(goto 3)
    # 3: first line of repeat info
    #   expect: repeat_info(store begin, goto 4) or ">Repeat"(goto 5)
    # 4: continuing to read repeat info
    #   expect: repeat_info(goto 4) or ">Repeat"(goto 5)
    # 5: part of MSA sequence
    #   expect: sequence(goto 6)
    # 6: After first sequence
    # expect: ">Repeat"(goto 5) or "REPEAT_TYPE"(return repeat, goto 2) or
    # "//"(return repeat, goto 0)

    state = 0
    region = RepeatRegion()
    identifier = ""

    def strip_comments(line):
        pat_comment = re.compile(r"\s*#.*$")
        return pat_comment.sub("", line)

    for i, line in enumerate(infile):
        line = strip_comments(line)
        if not line or line == '\n':
            continue

        LOG.debug("Line %d: %s", i, line[0:-1])

        if 0 == state:
            match = pat_identifier.match(line)
            if match:
                LOG.debug(
                    " *(0->1) Found identifier (storing \"%s\")",
                    match.group(1))

                identifier = match.group(1)
                state = 1
                continue

        elif 1 == state:
            match = pat_repeat_type.match(line)
            if match:
                LOG.debug(" *(1->2) Found REPEAT_TYPE")
                state = 2
                continue

            match = pat_protein_end.match(line)
            if match:
                LOG.debug(" *(1->0) Found protein end")
                state = 0
                continue

        elif 2 == state:
            # Entry Action: clear msa cache
            region = RepeatRegion(protein_id=identifier)

            match = pat_repeat_length.match(line)
            if match:
                LOG.debug(" *(2->3) Found REPEAT_LENGTH")
                state = 3
                continue

        elif 3 == state:
            match = pat_repeat_info.match(line)
            if match:
                LOG.debug(
                    " *(3->4) Found repeat_info (storing begin: \"%s\")",
                    match.group(1)
                )

                region.begin = int(match.group(1))
                state = 4
                continue

            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" *(3->5) Found repeat_header")
                state = 5
                continue

        elif 4 == state:
            match = pat_repeat_info.match(line)
            if match:
                LOG.debug(" *(4->4) Found repeat_info")
                state = 4
                continue

            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" *(4->5) Found repeat_header")
                state = 5
                continue

        elif 5 == state:
            match = pat_repeat_sequence.match(line)
            if match:
                LOG.debug(
                    " *(5->6) Found sequence (storing \"%s\")",
                    match.group(1))

                region.msa.append(match.group(1))

                state = 6
                continue

        elif 6 == state:
            match = pat_repeat_header.match(line)
            if match:
                LOG.debug(" *(6->5) Found repeat_header")
                state = 5
                continue
            match = pat_repeat_type.match(line)
            if match:
                LOG.debug(" *(6->2) Found REPEAT_TYPE, yielding.")

                state = 2
                if len(region.msa) >= 2:
                    yield region
                continue
            match = pat_protein_end.match(line)
            if match:
                LOG.debug(" *(6->0) Found protein end, yielding.")
                state = 0
                if len(region.msa) >= 2:
                    yield region
                continue

        else:
            raise AssertionError("Huh? Unknown parser state " + str(state))

################################## TRF - Benson ##########################


def trf_get_repeats(infile):
    r"""Read repeats from a TRF txt.html file stream file stream successively.

    Read repeats from a TRF txt.html file stream file stream successively.
    Postcondition: infile points to EOF.

    TRF output file syntax::

        Sequence: ``identifier``
             Indices: ``begin``--``end``
             \d [a-zA-Z]+
        #
             begin (repeat)*
             1  (consensus)*
        #
          (( \d (repeat)*
             \d  (consensus)*
          )?
             \d (repeat)*
             1  (consensus)*
          )+
             \d [a-zA-Z]+
        #
            ``Statistics``

    Args:
        infile (file stream): File stream from TRF output txt.html.
    (generated via e.g. ./trf404.linux64.exe FASTAFILE 2 7 7 80 10 50 500 -d > /dev/null
    If the -h flag is set, no .txt.html output is produced)

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout TRF output syntax.
    .. todo:: Does not search for the sequence identifier at current!
    """

    # find the name of the sequence ## CURRENTLY NOT IMPLEMENTED
    #pat_identifier = re.compile(r"Sequence: (\S+)")

    # find coordinates of the repeat region in the protein
    pat_coordinates = re.compile(r"Indices: (\d+)--(\d+)")

    # find a part of a repeat unit and its coordinate
    pattern_seq = re.compile(r"\s+(\d+) ([ACGT\- ]+)")

    # find the final tag 'Statistics'
    pat_statistics = re.compile(r"Statistics")

    # Our possible parser states:
    #
    # state 0: searching for identifier -> 1  # not necessary when sequence identifier is known
    # state 1: searching for repeat region coordinates -> 2
    # state 2: searching for beginning of MSA & save sequence to tmpMSA-> 4
    # state 3: new repeat unit: save sequence to tmpMSA -> 6
    # state 4: new repeat unit: save sequence to tmp_consensus -> 5
    # state 5:
    #           if sequence: save sequence to tmpMSA -> 6
    #           if end: save tmpMSA to preMSA; save tmp_consensus to consensus; Yield repeat region -> 1
    # state 6: check: new repeat unit? save tmp_consensus
    #  1: use last tmpMSA entry for new tmpMSA;
    #     save tmpMSA to preMSA;
    #     save tmp_consensus to consensus;
    #     save sequence to new tmp_consensus -> 3
    #  0: save sequence to tmp_consensus -> 5

    state = 1
    # identifier = ""  # Currently not implemented.
    preMSA = []
    consensus = []
    tmp_consensus = None
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])

        # CURRENTLY NOT IMPLEMENTED
        # if state == 0: #searching for sequenceMSA identifier
        #  tmp = pat_identifier.search(line)
        #  if tmp:
        #    state = 1
        #    identifier = tmp.group()
        #  continue
        # elif state == 1: #searching for TR boundaries (indices)

        if 1 == state:  # searching for repeat region coordinates
            search = pat_coordinates.search(line)
            if search:
                LOG.debug(" * (1->2) Found coordinates")
                state = 2
                region = RepeatRegion()
                region.begin = int(search.group(1))
                #region_end = search.group(2)
                short = False

        # searching for beginning of MSA & save sequence to tmpMSA-> 4
        elif state == 2:
            match = pattern_seq.match(line)
            if match and match.group(1) == str(region.begin):
                LOG.debug(" *(2->4) Found first row of first MSA repeat unit")
                state = 4
                if len(match.group(2).strip().split(" ")) > 1:
                    short = True
                    preMSA = match.group(2).strip().split(" ")
                    LOG.debug("Repeat unit is short; preMSA: %s", str(preMSA))
                else:
                    tmpMSA = [match.group(2).strip().split(" ")[0]]
                    LOG.debug(" tmpMSA: %s", str(tmpMSA))

        elif state == 3:  # new repeat unit: save sequence to tmpMSA -> 4
            match = pattern_seq.match(line)
            if match:
                LOG.debug(" *(3->5) Found first row of new repeat unit")
                state = 6
                if short:
                    preMSA += match.group(2).strip().split(" ")
                    LOG.debug("Repeat unit is short; preMSA: %s", str(preMSA))
                else:
                    tmpMSA.append(match.group(2).strip().split(" ")[0])
                    LOG.debug(" tmpMSA: %s", str(tmpMSA))
            # if end: save tmpMSA to preMSA; save tmp_consensus to consensus;
            # Yield repeat region -> 1
            if pat_statistics.search(line):
                LOG.debug(
                    " *(5->1) Encountered 'Statistics': No more repeats, yielding.")
                state = 1
                if not short:
                    preMSA.append("".join(tmpMSA))
                    consensus.append("".join(tmp_consensus))
                region.msa = getMSA(preMSA, consensus)
                if len(region.msa) >= 2:
                    yield region
                preMSA = []
                consensus = []

        elif state == 4:  # new repeat unit: save sequence to tmp_consensus -> 5
            match = pattern_seq.match(line)
            if match:
                LOG.debug(
                    " *(4->5) Found first consensus row of the repeat unit")
                state = 5
                if short:
                    consensus = match.group(2).strip().split(" ")
                    LOG.debug(
                        "Repeat unit is short;  consensus: %s",
                        str(consensus))
                else:
                    tmp_consensus = [match.group(2).strip().split(" ")[0]]
                    LOG.debug(" tmp_consensus: %s", str(tmp_consensus))

        elif state == 5:  # SEARCHING FOR MSA ROW
            # if sequence: save sequence to tmpMSA -> 6
            match = pattern_seq.match(line)
            if match:
                LOG.debug(" *(5->6) Found a MSA row")
                state = 6
                if short:
                    preMSA += match.group(2).strip().split(" ")
                    LOG.debug("Repeat unit is short; preMSA: %s", str(preMSA))
                else:
                    tmpMSA.append(match.group(2).strip().split(" ")[0])
                    LOG.debug(" tmpMSA: %s", str(tmpMSA))

            # if end: save tmpMSA to preMSA; save tmp_consensus to consensus;
            # Yield repeat region -> 1
            if pat_statistics.search(line):
                LOG.debug(
                    " *(5->1) Encountered 'Statistics': No more repeats, yielding.")
                state = 1
                if not short:
                    preMSA.append("".join(tmpMSA))
                    consensus.append("".join(tmp_consensus))
                region.msa = getMSA(preMSA, consensus)
                if len(region.msa) >= 2:
                    yield region
                preMSA = []
                consensus = []

        elif state == 6:  # new repeat unit? ## SEARCHING FOR CONSENSUS ROW
            match = pattern_seq.match(line)
            # 1: save tmp_consensus -> 3
            if match and match.group(1) == '1':
                LOG.debug(
                    " *(6->3) Found a consensus row of a new repeat unit")
                state = 3
                if short:  # NEEDS TO BE CODED
                    consensus += match.group(2).strip().split(" ")
                    LOG.debug(
                        "Repeat unit is short; consensus: %s",
                        str(consensus))
                else:
                    newMSA = tmpMSA.pop()
                    preMSA.append("".join(tmpMSA))
                    consensus.append("".join(tmp_consensus))
                    tmpMSA = [newMSA]
                    tmp_consensus = [match.group(2).strip().split(" ")[0]]

            # 0: save sequence to tmp_consensus -> 5
            elif match:
                LOG.debug(" *(6->5) Found a consensus row")
                state = 5
                tmp_consensus.append(match.group(2).strip().split(" ")[0])

            # YIELD
            # aha! there should have been a consensus sequence, but there is
            # not. Hence we are finished with this repeat!
            else:
                LOG.debug(' *(6->1) No consensus row: repeat finished')
                state = 1
                if short:
                    preMSA = preMSA[:-1]
                else:
                    preMSA.append(
                        "".join(
                            tmpMSA[
                                0:-
                                1]))  # The last tmpMSA entry was not a repeat unit
                    consensus.append("".join(tmp_consensus))
                LOG.debug(" preMSA: %s", str(preMSA))
                LOG.debug(" consensus: %s", str(consensus))
                region.msa = getMSA(preMSA, consensus)
                if len(region.msa) >= 2:
                    yield region
                preMSA = []
                consensus = []


def getMSA(sequenceMSA, consensusMSA):
    """ Derive the MSA from a strange combination of consensusMSA and sequenceMSA in TRF
    (Benson) txt.html output files

    Args:
        sequenceMSA (?):
        consensusMSA (?):

    Returns:
         msa (list of str): The multiple sequence alignment predicted by TRF.
    """

    msa = [""] * len(sequenceMSA)

    while consensusMSA:
        # CHECK for insertions
        insertion = 1
        while insertion and consensusMSA:
            insertion = 0
            for i_con in consensusMSA:
                if i_con and i_con[0] == "-":
                    insertion = 1
                    break
            # INCLUDE insertions into the msa
            if insertion:
                for i in range(len(consensusMSA)):
                    if consensusMSA[i] and consensusMSA[i][0] == "-":
                        msa[i] += sequenceMSA[i][0]
                        sequenceMSA[i] = sequenceMSA[i][1:]
                        consensusMSA[i] = consensusMSA[i][1:]
                    else:
                        msa[i] += "-"

        # CHECK for deletions and normal sequence
        if not consensusMSA[0]:
            break

        for i in range(len(consensusMSA)):
            # The last repeat unit can be shorter than the ones before
            if not sequenceMSA[i]:
                break

            msa[i] += sequenceMSA[i][0]
            sequenceMSA[i] = sequenceMSA[i][1:]
            consensusMSA[i] = consensusMSA[i][1:]

    return msa

# ############## HHrepID - Soeding ###########################################


def hhpredid_get_repeats(infile):
    r"""Read repeats from a HHREPID standard output (stdout) file stream successively.

    Read repeats from a HHREPID standard output (stdout) file stream successively.
    Postcondition: infile points to EOF.

    Layout of HHREPID standard output::

        protein ::=
             begin"-"\d    "+"\d repeatUnit
           ( \d"-"\d    "+"\d repeatUnit )+

    Args:
        infile (file stream): File stream from HHREPID standard output.
        [Generated by e.g.: ./hhrepid_32 -i FASTAFILE -v 0 -d cal.hhm -o INFILE]

    Returns:
         (Repeat): A generator function is returned that, when called in a loop,
         yields one repeat per iteration.

    .. todo:: Layout HHREPID output syntax.
    """

    # find a part of a repeat unit and its first coordinate
    # minus or \minus?

    pattern_repeat_unit_count = re.compile(r"Repeats\s+(\d+)")
    pattern_seq = re.compile(r"[A-Z]+(\d+).*(\d+)\-.*\+[\d]+ ([\-a-zA-Z.]+)")

    # Our possible parser states:

    # state1: Find number of repeat units n
    # state2: Find first (partial) row of the MSA
    # state3: Find all other (partial) rows of the MSA

    region = None
    state = 1
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])

        if 1 == state:  # Find 'Repeats' marker of new repeat
            search = pattern_repeat_unit_count.search(line)
            if search:
                LOG.debug(" *(1->2) Found repeat")
                state = 2
                n = int(search.group(1))

        elif 2 == state:  # Find first (partial) row of the MSA
            search = pattern_seq.search(line)
            if search:
                LOG.debug(" *(2->3) Found first repeat unit (part)")
                state = 3
                region = RepeatRegion()
                region.begin = int(search.group(2))
                region.msa = [""] * n
                region.msa[int(search.group(1)) -
                           1] = search.group(3).replace('.', '-').upper()

        elif 3 == state:  # Find all other (partial) rows of the MSA
            search = pattern_seq.search(line)
            if search:
                LOG.debug(" *(3->3) Found other repeat unit (part)")
                region.msa[int(search.group(1)) -
                           1] += search.group(3).replace('.', '-').upper()
            else:
                search = pattern_repeat_unit_count.search(line)
                if search:
                    LOG.debug(" *(3->2) Yield Repeat, begin next")
                    state = 2
                    n = int(search.group(1))
                    if len(region.msa) >= 2:
                        yield region
                        region = None
                    else:
                        LOG.warning(
                            "HHPREDID: Msa too short %s", str(
                                region.msa))

    # Yield final repeat region.
    if region is not None:
        if len(region.msa) >= 2:
            yield region
        else:
            LOG.warning("HHPREDID: Msa too short %s", str(region.msa))

####################################### Phobos TRF  ######################


# def phobos_get_repeats(infile):
#     """ Read repeats from a PHOBOS output file stream successively.

#     Read repeats from a PHOBOS output file stream successively.
#     Postcondition: infile points to EOF.

#     Args:
#         infile (file stream): File stream from PHOBOS output.

#     Returns:
#          (Repeat): A generator function is returned that, when called in a loop,
#          yields one repeat per iteration.

#     .. todo:: Show PHOBOS output syntax.
#     """

#     pattern_begin = re.compile(r"(\d+) :\s+\d")
#     pattern_seq = re.compile(r"([\-ACGT]+)")

#     # Our possible parser states:
#     #
#     # state 1: Find TR begin
#     # state 2: Find first repeat unit
#     # state 3: Find repeat units

#     state = 1
#     for i, line in enumerate(infile):
#         LOG.debug("Line %d: %s", i, line[0:-1])
#         if 1 == state:  # Find TR offset
#             search = pattern_begin.search(line)
#             if search and search.groups()[0] is not None:
#                 LOG.debug(" *(1->2) Found tandem repeat begin")
#                 state = 2
#                 region = RepeatRegion()
#                 region.begin = int(search.groups()[0])
#                 region.msa = []

#         elif 2 == state:  # Find all other repeat units
#             match = pattern_seq.search(line)
#             if match and match.groups()[0] is not None:
#                 LOG.debug(" *(2->3) Found first repeat unit")
#                 region.msa.append(match.groups()[0])
#                 state = 3

#         elif 3 == state:  # Find all other repeat units
#             match = pattern_seq.search(line)
#             if match and match.groups()[0] is not None:
#                 LOG.debug(" *(3->3) Found a repeat unit")
#                 region.msa.append(match.groups()[0])
#             else:
#                 LOG.debug(" *(3->1) repeat region finished, yielding.")
#                 state = 1
#                 if len(region.msa) >= 2:
#                     yield region
#                 else:
#                     LOG.warning("phobos: Msa too short %s", str(region.msa))

def phobos_get_repeats(self, infile):
    """ Read repeats from a PHOBOS output file stream successively.
    Read repeats from a PHOBOS output file stream successively.
    Postcondition: infile points to EOF.
    Args:
        infile (file stream): File stream from PHOBOS output.
    Returns:
        (Repeat): A generator function is returned that, when called in a loop,
        yields one repeat per iteration.
    .. todo:: Show PHOBOS output syntax.
    """

    pattern_begin = re.compile(r"(\d+) :\s+\d")
    pattern_seq = re.compile(r"(?<![0-9])([\-ACGT]+)")
    pattern_actual_tr = re.compile(r"(?<![0-9])([\-ACGTN]+)")
    pattern_indels = re.compile(r"(?:\(([0-9]+,[DI])\))")

    # Our possible parser states:
    #
    # state 1: Find TR begin in sequence
    # state 2: Find insertions
    # state 3: Find actual TR sequence
    # state 4: Find ideal TR sequence, extract repeat units and yield

    state = 1
    for i, line in enumerate(infile):
        LOG.debug("Line %d: %s", i, line[0:-1])
        if 1 == state:  # Find TR offset
            search = pattern_begin.search(line)
            if search and search.groups()[0] is not None:
                LOG.debug(" *(1->2) Found tandem repeat begin")     
                # Find TR begin in sequence, init RepeatRegion               
                region = RepeatRegion()
                region.begin = int(search.groups()[0])
                region.msa = []
                # Find TR concensus unit
                unit_match = pattern_seq.search(line)
                unit = unit_match.groups()[0]
                # Update state
                state = 2

        elif 2 == state:  # Find potential insertions
            # Extract region insertions and deletions as list of tuples ([0-9]+, [DI])
            indels = [(int(i.split(",")[0]), i.split(",")[1]) for i in pattern_indels.findall(line)]
            del_counter = 0
            insertion_sites = []
            # The position of an insertion needs to be incremented once for each deletion that
            ## occurs to the left of it in the repeat region
            for indel in indels:
                if indel[1] == "D":
                    del_counter += 1
                    # continue
                elif indel[1] == "I":
                    insertion_sites.append(int(indel[0]) - 1 + del_counter)
            state = 3

        elif 3 == state:  # Find actual TR sequence
            # match = pattern_seq.search(line)
            match = pattern_actual_tr.search(line)
            if match and match.groups()[0] is not None:
                LOG.debug(" *(2->3) Found actual repeat region")
                actual_tr = match.groups()[0]
                state = 4

        elif 4 == state:  # Find ideal TR sequence, extract repeat units and yield
            match = pattern_seq.search(line)
            if match and match.groups()[0] is not None:
                LOG.debug(" *(3->3) Found ideal repeat region")
                ideal_tr = match.groups()[0]
                
                # get locations of unit start and end sites, as well as gap position relative to unit starts
                unit_locations_raw, gaps = self.map_unit_boundaries(ideal_tr, unit, insertion_sites)

                # remove out of range units, trim border-overlapping units
                trimmed_unit_locations = self.trim_unit_locations(unit_locations_raw, len(unit), len(actual_tr))                    
                # extract units from TR to construct msa      
                for i in trimmed_unit_locations:
                    trimmed_unit_locations[i]["unit"] = actual_tr[trimmed_unit_locations[i]["begin"]:trimmed_unit_locations[i]["end"]+1]                    
                
                # Get first unit, fill out by adding "-" to start of unit until it is as long as the consensus unit
                lead_unit = list(trimmed_unit_locations.items())[0][1]
                if len(lead_unit["unit"]) != len(unit):
                    difference = len(unit) - len(lead_unit["unit"]) 
                    lead_unit["unit"] = "-" * difference + lead_unit["unit"]
                # Get last unit, fill out by adding "-" to end of unit until it is as long as the consensus unit
                tail_unit = list(trimmed_unit_locations.items())[-1][1]
                if len(tail_unit["unit"]) != len(unit):
                    difference = len(unit) - len(tail_unit["unit"]) 
                    tail_unit["unit"] = tail_unit["unit"] + "-" * difference 

                # insert "-" at the gap locations in each unit
                ## make sure to not add "-" if the insert causing the gap is in the current unit
                if gaps:
                    for i in trimmed_unit_locations: # will return keys
                        counter = 0
                        # for gap in unit_gaps:
                        for gap in gaps:
                            if gap in trimmed_unit_locations[i]["inserts"]:
                                # unit contains gap-causing insert, do not add "-" to this unit
                                counter += 1
                                continue
                            # correct for potential insertion of earlier gap in unit
                            gap += counter
                            counter += 1
                            # insert "-" at proper position in unit
                            trimmed_unit_locations[i]["unit"] = trimmed_unit_locations[i]["unit"][:gap] + \
                                "-" + trimmed_unit_locations[i]["unit"][gap:]
                # generate msa list from OrderedDict, msa_list is list of strings representing the repeat units
                msa_list = [trimmed_unit_locations[i]["unit"] for i in trimmed_unit_locations]

                LOG.debug(" *(3->1) repeat region finished, yielding.")
                state = 1                
                if len(msa_list) >= 2:
                    region.msa = msa_list
                    yield region
                else:                        
                    LOG.warning("phobos: Msa too short %s", str(region.msa))

def map_unit_boundaries(self, ideal_tr, unit, inserts):   
    """ Map rough unit begin and start sites on repeat region
    Start by finding a match for the consensus TR unit in ideal TR representation provided by PHOBOS.
    Starting from this unit, generate a left and a right flank of units until whole repeat region is covered. Merge 
    these three components and return as one OrderedDict
    This will likely contain units that extend beyond the region begin and end, which will be trimmed later
    using trim_unit_locations(). Unit begin and start sites for left and right flank will be corrected for
    insertions using correct_inserts_left() and correct_inserts_right(), respectively.

    Parameters
    ideal_tr (str):     String representation of the ideal repeat region, extracted from PHOBOS output
    unit (str):         String representation of consensus TR unit, extracted from PHOBOS output
    inserts (list(int)):     
                        List of integers specifying where in the repeat region insertions have occurred,
                        extracted from PHOBOS output (preprocessed to account for deletions prior to inserts)

    Returns
    merged_odict (OrderedDict):
                        An ordered dictionary where the keys are (ascending) integers, and values are dictionaries
                        describing the begin, end and insertion sites for each unit
    gaps (list(int)):   A list of gap positions (relative to unit start site) that will be needed to insert gaps
                        in proper places later
    """     
    region_l = len(ideal_tr)
    unit_l = len(unit)
    i = 0
    first_match = ""
    # Slide consensus unit over ideal sequence until match is found. Extract the boundaries of this match
    ## and construct 'first_match' tuple
    while i <= region_l:
        if ideal_tr[i:i + unit_l] == unit:
            first_match = (i, i + unit_l - 1)
            break
        i += 1
    if not first_match:
        raise ValueError("No match for unit '{}' found in ideal TR region".format(unit))                

    # initialize and populate left flank as (begin, end)
    ## keep adding units to the left of first_match until out of range of region
    left_flank = [(first_match[0] - unit_l, first_match[1] - unit_l)]
    while left_flank[0][1] >= 0:
        new_unit =(
            left_flank[0][0] - unit_l,
            left_flank[0][1] - unit_l
        )
        # insert unit tuples at begin of list so they remain in order of appearance in repeat region
        left_flank.insert(0, new_unit)

    # initialize and populate right flank
    ## keep adding units to the right of first_match until out of range of region
    right_flank = [(first_match[0] + unit_l, first_match[1] + unit_l)]
    while right_flank[-1][0] <= region_l - 1:
                new_unit = (
                    right_flank[-1][0] + unit_l, 
                    right_flank[-1][1] + unit_l
                )
                right_flank.append(new_unit)

    # Now that all (begin, end) tuples have been generated, an OrderedDict is made for both flanks and for the first match.
    ## Keys will be ascending integers and values will be dicts {begin: int, end: int, inserts: list}
    key_counter = 1
    left_flank_odict = OrderedDict()
    for i in left_flank:
        left_flank_odict[key_counter] = {
            "begin": i[0],
            "end": i[1],
            "inserts": []
        }
        key_counter += 1

    first_match_odict = OrderedDict([(key_counter, {"begin": first_match[0], "end": first_match[1], "inserts": []})])
    key_counter += 1

    right_flank_odict = OrderedDict()
    for i in right_flank:
        right_flank_odict[key_counter] = {
            "begin": i[0],
            "end": i[1],
            "inserts": []
        }
        key_counter += 1
    
    # Modify unit begin and end sites to correct for insertions
    ## Also extract gap postions relative to unit begin for later use
    all_gaps = None
    if inserts:
        left_flank_odict, left_gaps = self.correct_inserts_left(left_flank_odict, inserts)
        right_flank_odict, right_gaps = self.correct_inserts_right(right_flank_odict, inserts)
        all_gaps = left_gaps.union(right_gaps)
        all_gaps = sorted(list(all_gaps))

    # Merge the three OrderedDicts and return
    ## Note: out of range units are still included here, they will be removed/trimmed later
    merged_odict = OrderedDict(list(left_flank_odict.items()) + list(first_match_odict.items()) + list(right_flank_odict.items()))
    return merged_odict, all_gaps

def correct_inserts_left(self, left_flank, insertion_sites):  
    """ Correct for insertion sites occurring in the left flank of TR region.
    When an insertion occurs, the boundaries of units in the left flank have to be corrected and shifted
    toward the begin of the repeat region

    Parameters
    left_flank (OrderedDict):      
                            Ordered dictionary with ascending integers as keys and dictionaries describing repeat units as values.
                            e.g. {1:{begin: int, end: int, insertions: list(int)}} Unit begin and end positions are to be
                            corrected for insertions
    insertion_sites (list):
                            List of integers specifying where in the repeat region insertions have occurred

    Returns
    left_flank (OrderedDict):
                            Same Ordered dictionary as input, but with the corrected begin and end position of TR units
    gap_positions (set):    A set of relative gap positions. Relative meaning that gap position are counted from unit begin site,
                            not region begin site
    """      
    gap_positions = set()
    for insert in insertion_sites:
        if insert > list(left_flank.items())[-1][1]["end"]: # insert is in right flank, skip                
            continue
        for i in range(list(left_flank.items())[-1][0], 0, -1):
            if insert >= left_flank[i]["begin"] and insert <= left_flank[i]["end"]:
                # Insert in unit: push begin back one position, update gaps                    
                left_flank[i]["begin"] -= 1                       
                insert_site = insert - left_flank[i]["begin"]                                     
                left_flank[i]["inserts"].append(insert_site)
                gap_positions.add(insert_site)                    
            elif insert > left_flank[i]["end"]:
                # Insert ahead of unit, push both begin and end back 
                left_flank[i]["begin"] -= 1
                left_flank[i]["end"] -= 1      
    return left_flank, gap_positions

def correct_inserts_right(self, right_flank, insertion_sites):
    """ Correct for insertion sites occurring in the right flank of TR region.
    When an insertion occurs, the boundaries of units in the right flank have to be corrected and shifted
    toward the end of the repeat region

    Parameters
    right_flank (OrderedDict):      
                            Ordered dictionary with ascending integers as keys and dictionaries describing repeat units as values.
                            e.g. {1:{begin: int, end: int, insertions: list(int)}} Unit begin and end positions are to be
                            corrected for insertions
    insertion_sites (list):
                            List of integers specifying where in the repeat region insertions have occurred

    Returns
    right_flank (OrderedDict):
                            Same Ordered dictionary as input, but with the corrected begin and end position of TR units
    gap_positions (set):    A set of relative gap positions. Relative meaning that gap position are counted from unit begin site,
                            not region begin site
    """
    gap_positions = set()
    for insert in insertion_sites:
        if insert < list(right_flank.items())[0][1]["begin"]: # insert is in left flank, skip                
            continue
        for i in range(list(right_flank.items())[0][0], list(right_flank.items())[-1][0] + 1):
            if insert >= right_flank[i]["begin"] and insert <= right_flank[i]["end"]:
                # Insert in unit: push end up one position, update gaps
                right_flank[i]["end"] += 1                        
                insert_site = insert - right_flank[i]["begin"]                                    
                right_flank[i]["inserts"].append(insert_site)
                gap_positions.add(insert_site)                    
            elif insert < right_flank[i]["begin"]:
                # Insert prior to unit, push both begin and end up
                right_flank[i]["begin"] += 1
                right_flank[i]["end"] += 1             
    return right_flank, gap_positions

def trim_unit_locations(self, units_odict, unit_length, region_length):
    """ Correct units that overlap repeat region begin or end (e.g. (-3, 4) -> (0, 4)) 
    and trim units that fall out of the repeat region range (e.g. remove (-10, -3))

    Parameters
    units_odict (OrderedDict):       
                            Ordered dictionary with ascending integers as keys and dictionaries describing repeat units as values.
                            e.g. {1:{begin: int, end: int, insertions: list(int)}}
                            Note: unit boundaries should have been corrected for insertions by the
                            appropriate functions previously
    unit_length (int):      Length of the consensus repeat unit
    region_length (int):    Length of the repeat region

    Returns
    units_odict (OrderedDict):
                            Processed input OrderedDict: Out of range units have been removed and 
                            region boundary overlaps have been corrected
    """
    for i in range(1, len(units_odict) + 1):
        if units_odict[i]["end"] < 0:
            del units_odict[i]
        elif units_odict[i]["begin"] >= region_length:
            del units_odict[i]
        elif units_odict[i]["begin"] < 0:
            units_odict[i]["begin"] = 0
        elif units_odict[i]["end"] > region_length:
            units_odict[i]["end"] = region_length
    return units_odict
