import numpy

class MdpaWriter():
    def __init__(self, name, data):
        self.fid = open(name + ".mdpa", 'w')
        self.data = data
        self.fid.write("//KRATOS isogeometric application data file\n")
        self.fid.write("//(c) 2014 Hoang Giang Bui, Ruhr-University Bochum\n\n")

    def Finalize(self):
        self.fid.close()

    def WriteMdpaData(self):
        self.fid.write("Begin ModelPartData\n")
        self.fid.write("End ModelPartData\n\n")

    def WriteProperties(self):
        cnt = 1
        for str_layer in self.data.layer_list:
            self.data.layer_attributes[str_layer]['id'] = cnt
            self.fid.write("Begin Properties " + str(cnt) + "\n")
            for key in self.data.layer_attributes[str_layer]['properties']:
                self.fid.write(str(key) + " " + str(self.data.layer_attributes[str_layer]['properties'][key]) + "\n")
            self.fid.write("End Properties\n\n")
            cnt = cnt + 1

    def WriteNodes(self):
        self.fid.write("Begin Nodes\n")
        for str_layer in self.data.layer_list:
            for i_node in self.data.layer_nodes_sets[str_layer]:
                n = self.data.layer_nodes_sets[str_layer][i_node]
                self.fid.write(str(i_node) + " " + str(n[0]) + " " + str(n[1]) + " " + str(n[2]) + "\n")
        self.fid.write("End Nodes\n\n")

    def WriteEntities(self):
        for str_layer in self.data.layer_list:
            if self.data.layer_attributes[str_layer]['type'] == "element" or \
                    self.data.layer_attributes[str_layer]['type'] == "Element" or \
                    self.data.layer_attributes[str_layer]['type'] == "Elements":
                str_entity = "Elements"
                str_entity_data = "ElementalData"
            elif self.data.layer_attributes[str_layer]['type'] == "condition" or \
                    self.data.layer_attributes[str_layer]['type'] == "Condition" or \
                    self.data.layer_attributes[str_layer]['type'] == "Conditions":
                str_entity = "Conditions"
                str_entity_data = "ConditionalData"
            else:
                print("Error: unknown entity type")
                exit(0)

            # write connectivities
            self.fid.write("Begin " + str_entity + " " + str(self.data.layer_attributes[str_layer]['name']) + "\n")
            for i_entity in self.data.layer_entities_sets[str_layer]:
                self.fid.write(str(i_entity) + " " + str(self.data.layer_attributes[str_layer]['id']))
                for i_node in self.data.layer_entities_sets[str_layer][i_entity]:
                    self.fid.write(" " + str(i_node))
                self.fid.write("\n")
            self.fid.write("End " + str_entity + "\n\n")

            # write entity data
            for str_info in self.data.layer_entity_info_sets[str_layer]:
                self.fid.write("Begin " + str_entity_data + " " + str(str_info) + "\n")
                for i_entity in self.data.layer_entity_info_sets[str_layer][str_info]:
                    temp = self.data.layer_entity_info_sets[str_layer][str_info][i_entity]
                    if isinstance(temp, list):
                        temp = numpy.array(temp)
                        dim = len(temp.shape)
                        if dim == 1:
                            self.fid.write(str(i_entity) + " " + str(self.GetVectorString(temp)) + "\n")
                        elif dim == 2:
                            self.fid.write(str(i_entity) + " " + str(self.GetMatrixString(temp)) + "\n")
                        else:
                            pass
                    else:
                        self.fid.write(str(i_entity) + " " + str(temp) + "\n")
                self.fid.write("End " + str_entity_data + "\n\n")

    def GetVectorString(self, v):
        st = '[' + str(len(v)) + '] (';
        for i in range(0, len(v)):
            if i != len(v)-1:
                st = st + str(v[i]) + ", "
            else:
                st = st + str(v[i]) + ")"
        return st

    def GetMatrixString(self, v):
        st = '[' + str(v.shape[0]) + ', ' + str(v.shape[1]) + '] ((';
        for i in range(0, len(v[0])):
            if i != len(v[0])-1:
                st = st + str(v[0][i]) + ", "
            else:
                st = st + str(v[0][i]) + "),("
        for i in range(0, len(v[1])):
            if i != len(v[1])-1:
                st = st + str(v[1][i]) + ", "
            else:
                st = st + str(v[1][i]) + "))"
        return st

