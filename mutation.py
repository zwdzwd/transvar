import re, sys

def parse_mutation_str(mut_str):

    mp = re.match(r'(p.)?([A-Z*]?)(\d+)([A-Z*]?)', mut_str)
    mn = re.match(r'(c.)?(\d+)([ATGC]?)>([ATGC]?)', mut_str)

    if mp and not mn:
        is_codon = True
        ref = mp.group(2) if mp.group(2) else '.'
        pos = int(mp.group(3))
        alt = mp.group(4) if mp.group(4) else '.'
    elif mn and not mp:
        is_codon = False
        pos = int(mn.group(2))
        ref = mn.group(3) if mn.group(3) else '.'
        alt = mn.group(4) if mn.group(4) else '.'
    else:
        sys.stderr.write('Cannot infer mutation type for %s' % mut_str)

    return (is_codon, pos, ref, alt)
