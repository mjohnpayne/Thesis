__author__ = 'mjohnpayne'

import sys

infile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/node_assignments_tree.nwk','r')
outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/node_assignments_tree_nos_only.nwk','w')
infile = infile.read()

print infile

def rem_names(s):
    if "_" in s:
        p = s.find("_")
        nstr = s[:p-3]+s[p+1:]
        return rem_names(nstr)
    else:
        return s

def rem_supports(s):
    if ":" in s:
        p = s.find(":")
        nstr = s[:p]+s[p+11:]
        return rem_supports(nstr)
    else:
        return s

def swap_no_for_name(s):
    if "_" in s:
        p = s.find("_")
        nstr = s[p-3:]
        n = nstr.find(':')
        id = nstr[:n]
        id = id[4:] + ':' + id[:3]
        print id
        nstr = s[:p-3]+s[p+11:]
        return rem_supports(nstr)
    else:
        return s

swap_no_for_name(infile)


out1 = rem_names(infile)

print out1
#
# out2 = rem_supports(out1)
#
# print out2
#
# out3 = rem_supports(infile)
#
# print out3
#
outfile.write(out1)
#
# outfile.close()