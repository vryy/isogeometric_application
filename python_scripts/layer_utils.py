import math
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.DiscontinuitiesApplication import *
from KratosMultiphysics.IsogeometricApplication import *

#
# Collapse each layer in Layers; every layer maintains its local numbering and connectivities
#
def CollapseLocal(Layers):
    tol = 1.0e-6

    # collapse each layer specified in layer list
    for str_layer in Layers.layer_list:
        node_map = {}
        non_repetitive_nodes = []
        # iterate all the nodes to check for node repetition
        for i_node in self.layer_nodes_sets[str_layer]:
            n = self.layer_nodes_sets[str_layer][i_node]
            for i_other_node in self.layer_nodes_sets[str_layer]:
                n1 = self.layer_nodes_sets[str_layer][i_other_node]
                d = math.sqrt(math.pow(n[0] - n1[0], 2) + math.pow(n[1] - n1[1], 2) + math.pow(n[2] - n1[2], 2))
                if d < tol:
                    node_map[i_node] = i_other_node
                    if i_other_node == i_node:
                        non_repetitive_nodes.append(i_node)
                    break

        # reform the layer nodal set
        new_nodes_set = {}
        node_map_nonrepetitive = {}
        cnt = 1
        for i_node in non_repetitive_nodes:
            new_nodes_set[cnt] = self.layer_nodes_sets[str_layer][i_node]
            node_map_nonrepetitive[i_node] = cnt
            cnt = cnt + 1
        self.layer_nodes_sets[str_layer] = new_nodes_set

        # reform the entity connectivities
        for i_entity in self.layer_entities_sets[str_layer]:
            new_entity = []
            for i_node in self.layer_entities_sets[str_layer][i_entity]:
                new_entity.append(node_map_nonrepetitive[node_map[i_node]])
            self.layer_entities_sets[str_layer][i_entity] = new_entity

    # turn on local collapse flag
    Layers.is_collapse_local = True

#
# Collapse all layers in Layers; all the layers nodal set will be reformed after collapsing. For this to work correctly, the model needs to be collapsed locally first to keep a list of non-repeated nodes
#
def Collapse(Layers):
    if not Layers.is_collapse_local == False:
        print("Error: the layers are not locally collapsed yet")
        sys.exit(0)

    tol = 1.0e-6
    binning_util = SpatialBinningUtility(0, 0, 0, 0.5, 0.5, 0.5, tol)

    # extract all the nodes and put into the spatial binning
    node_map = {}
    for str_layer in Layer.layer_list:
        node_map[str_layer] = {}
        for i_node in Layer.layer_nodes_sets[str_layer]:
            n = Layer.layer_nodes_sets[str_layer][i_node]
            new_id = binning_util.AddNode(n[0], n[1], n[2])
            node_map[str_layer][i_node] = new_id

