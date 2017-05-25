__author__ = 'mjohnpayne'

from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace
import math


#
def my_layout(node):
    if node.name != "45":
        if node.is_leaf():
             # If terminal node, draws its name
            name_face = AttrFace("species", fsize=16)
            pie = PieChartFace([changes[node.name][0],changes[node.name][1]],changes[node.name][2],changes[node.name][2],["Green","Red"])
            pie.opacity = 0.5
            #name_face = AttrFace("name", fsize=16)
            faces.add_face_to_node(name_face, node, column=0, position="branch-right")
            faces.add_face_to_node(pie, node, column=0, position="float")
        else:
             # If internal node, draws label with smaller font size
            #name_face = AttrFace("name", fsize=10)
            pie = PieChartFace([changes[node.name][0],changes[node.name][1]],changes[node.name][2],changes[node.name][2],["Green","Red"])
            pie.opacity = 0.5
            #faces.add_face_to_node(name_face, node, column=0, position="branch-right")
            faces.add_face_to_node(pie, node, column=0, position="float")


ts = TreeStyle()
# Do not add leaf names automatically ts.show_leaf_name = False
# Use my custom layout
ts.show_leaf_name = False
ts.layout_fn = my_layout





t = PhyloTree('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/node_assignments_tree_nos_only_dbl.nwk', format=1)

dataorder = ['A. fumigatus','N. fisheri','A. clavatus','A. terreus','A. flavus','A. oryzae','A. niger','A. nidulans','P. decumbens','P. roquefortii','P. chrysogenum','P. digitatum','T. stipitatus','P. funiculosum','T. marneffei','T. flavus','A. dermatiditis','H. capsulatum','P. brasiliensis','C. immitis','U. reesei','T. equinum','T. tonsurans']
nos = ["1",'2','4','6','7','8','12','14','16','17','18','19','24','25','26','27','32','33','35','37','38','40','41']

branch_to_node = {24:25,39:38,41:43,33:35,5:5,31:0,18:22,14:16,40:42,42:41,35:36,36:32,27:30,15:2,26:29,12:15,29:26,21:19,11:4,32:34,6:11,17:20,22:17,16:18,13:3,34:33,43:37,3:6,7:13,8:14,37:39,10:10,44:31,30:24,20:21,2:8,1:7,38:40,28:28,4:9,25:27,19:23,23:1,9:12}
#print branch_to_node

inchanges = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/gene_gain_loss_table_no_singles.txt','r')

changes = {}
for i in inchanges:
    if '#' not in i:
        inf = i[i.find(">")+1:].strip('\n').split('\t')
        gain = int(inf[1])
        loss = int(inf[2])
        tot = gain+loss
        pgain = gain/float(tot)*100
        ploss = 100-pgain
        size = tot
        size = math.sqrt((size/math.pi))
        changes[inf[0]] = (pgain,ploss,size)

print changes




dt = {}

for i in range(len(nos)):
    dt[nos[i]] = dataorder[i]

def get_species_name(node_name_string):
    # Species code is the first part of leaf name (separated by an # underscore character)
    spcode = node_name_string
    # We could even translate the code to complete names
    code2name = dt
    return code2name[spcode]

t.set_species_naming_function(get_species_name)

t.show(tree_style=ts)
#t.render("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/eurot_gain_loss_tree_nosingles_dbl.pdf",tree_style=ts)
