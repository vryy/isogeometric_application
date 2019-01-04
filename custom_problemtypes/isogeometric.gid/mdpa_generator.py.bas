import re
import math

def CheckSimplePattern(arg, pattern):
    collect = set()
    if(len(arg) > len(pattern)):
        for i in range(0, len(arg)):
            if(arg[i] in pattern):
                collect.add(arg[i])
        if(pattern == collect):
            params = {}
            for i in range(0, len(arg)):
                if(arg[i] in pattern):
                    params[arg[i]] = arg[i + 1]
            return [True, params]
        else:
            return [False]
    else:
        return [False]

def CheckLinePattern(arg):
    if(arg[0] == 'Line'):
        params = {}
        params[arg[0]] = arg[1]
        params[arg[2]] = arg[3]
        return [True, params]
    else:
        return [False]

def CheckPointPattern(arg):
    if(arg[0] == 'Point'):
        params = {}
        params['iu'] = arg[2]
        params['iv'] = arg[1]
        params[arg[3]] = []
        for i in range(4, len(arg)):
            params[arg[3]].append(arg[i])
        return [True, params]
    else:
        return [False]

def CheckKnotPattern(arg):
    pattern = {'Number', 'of', 'knots', 'in'}
    collect = set()
    if(len(arg) > len(pattern)):
        for i in range(0, len(arg)):
            if(arg[i] in pattern):
                collect.add(arg[i])
        if(pattern == collect):
            if(arg[4] == 'U'):
                return [True, 'U Knots', int(arg[5])]
            elif(arg[4] == 'V'):
                return [True, 'V Knots', int(arg[5])]
        else:
            return [False]
    else:
        return [False]

def CheckWeightPattern(arg):
    pattern = {'Rational', 'weights'}
    collect = set()
    if(len(arg) > len(pattern)):
        for i in range(0, len(arg)):
            if(arg[i] in pattern):
                collect.add(arg[i])
        if(pattern == collect):
            return [True]
        else:
            return [False]
    else:
        return [False]

def CheckDegreePattern(arg):
    pattern = {'Number', 'of', 'Control', 'Points', 'Degree'}
    collect = set()
    if(len(arg) > len(pattern)):
        for i in range(0, len(arg)):
            if(arg[i] in pattern):
                collect.add(arg[i])
        if(pattern == collect):
            params = {}
            params['U Dim'] = arg[4]
            params['V Dim'] = arg[5]
            params['U Degree'] = arg[7]
            params['V Degree'] = arg[8]
            return [True, params]
        else:
            return [False]
    else:
        return [False]

def ReadRawData(filename):
    fid = open(filename, 'r')

    NURBSURFACE_KW = 'NURBSURFACE'
    END_KW = 'END'

    nurbs_data = []

    lines = fid.readlines()
    line_read = 0

    while (line_read < len(lines) ):
        line = lines[line_read]
        tmp = line.split()
        
        if(len(tmp) == 0):
            line_read = line_read + 1
            continue
            
        if(tmp[0] == NURBSURFACE_KW):
            nurbs = {}
        
            #handle the nurbs surface data
            line_read = line_read + 1
            line = lines[line_read]
            tmp = re.split(r"[: =,]+", line.strip())
            while(True):
                if(len(tmp) == 0):
                    continue
                if(tmp[0] == END_KW):
                    break
                
                pattern = {'Num', 'HigherEntity', 'conditions', 'material'}
                res = CheckSimplePattern(tmp, pattern)
                if(res[0] == True):
                    for key in res[1]:
                        nurbs[key] = res[1][key]
                
                pattern = {'LAYER'}
                res = CheckSimplePattern(tmp, pattern)
                if(res[0] == True):
                    for key in res[1]:
                        nurbs[key] = res[1][key]
                
                pattern = {'NumLines'}
                res = CheckSimplePattern(tmp, pattern)
                if(res[0] == True):
                    for key in res[1]:
                        nurbs[key] = res[1][key]
                
                res = CheckLinePattern(tmp)
                if(res[0] == True):
                    if 'Lines' in nurbs:
                        nurbs['Lines'].append(res[1])
                    else:
                        nurbs['Lines'] = []
                        nurbs['Lines'].append(res[1])
                
                res = CheckPointPattern(tmp)
                if(res[0] == True):
                    if 'Points' in nurbs:
                        nurbs['Points'].append(res[1])
                    else:
                        nurbs['Points'] = []
                        nurbs['Points'].append(res[1])
                
                res = CheckKnotPattern(tmp)
                if(res[0] == True):
                    nurbs[res[1]] = []
                    for i in range(0, res[2]):
                        line_read = line_read + 1
                        line = lines[line_read]
                        tmp = re.split(r"[: =,]+", line.strip())
                        nurbs[res[1]].append(tmp[3])
                
                res = CheckWeightPattern(tmp)
                if(res[0] == True):
                    nurbs['Weights'] = []
                    for i in range(0, len(nurbs['Points'])):
                        line_read = line_read + 1
                        line = lines[line_read]
                        tmp = re.split(r"[: =,]+", line.strip())
                        nurbs['Weights'].append(tmp[0])
                
                res = CheckDegreePattern(tmp)
                if(res[0] == True):
                    for key in res[1]:
                        nurbs[key] = res[1][key]
                            
                line_read = line_read + 1
                line = lines[line_read]
                tmp = re.split(r"[: =,]+", line.strip())
            
            nurbs_data.append(nurbs)
        line_read = line_read + 1
    fid.close()
    model_data = {}
    model_data['nurbs_data'] = nurbs_data
    return model_data

def CheckCoincidentNode(node, nodes_list):
    for i in range(0, len(nodes_list)):
        distance = 0.0
        for j in range(0, len(node)):
            distance = distance + math.pow(float(node[j]) - float(nodes_list[i][j]), 2)
        distance = math.sqrt(distance)
        if distance < 1.0e-6:
            return i + 1
    return -1
    
def GenerateHeader(fid):
    fid.write('//KRATOS isogeometric application data file\n')
    fid.write("//(c) 2013 Hoang Giang Bui, Ruhr-University Bochum\n\n")

    fid.write('Begin ModelPartData\n')
    fid.write('//  VARIABLE_NAME value\n')
    fid.write('End ModelPartData\n\n')

def GenerateProperties(fid, model_data):
    props = set()
    for i in range(0, len(model_data['nurbs_data'])):
        props.add(int(model_data['nurbs_data'][i]['material']))
    
    #todo: extract data for material
    
    cnt = 1
    materials_map = {}
    for mat in props:
        if mat == 4: # UserDefined material
            fid.write('Begin Properties ' + str(cnt) + '\n')
            fid.write("GRAVITY [3] ( 0.0, 0.0, 0.0 )\n")
            fid.write('End Properties\n\n')
        elif mat == 2: # Isotropic3D material
            fid.write('Begin Properties ' + str(cnt) + '\n')
            fid.write('DENSITY 0.0\n')
            fid.write('BODY_FORCE [3] (0.0, 0.0, 0.0)\n')
            fid.write("GRAVITY [3] ( 0.0, 0.0, 0.0 )\n")
            fid.write('YOUNG_MODULUS 2.0e9\n')
            fid.write('POISSON_RATIO 0.3\n')
            fid.write('End Properties\n\n')
        materials_map[mat] = cnt
        cnt = cnt + 1
    model_data['materials_map'] = materials_map
    model_data['materials_list'] = props
            
def GenerateMdpa(filename, model_data):
    fid = open(filename, 'w')
    
    GenerateHeader(fid)
    GenerateProperties(fid, model_data)
    
    #extract all the control points within the model
    nodes = []
    nurbs_indices = []
    for i in range(0, len(model_data['nurbs_data'])):
        node_indices = []
        for node in model_data['nurbs_data'][i]['Points']:
            i = CheckCoincidentNode(node['coords'], nodes)
            if(i == -1):
                #new node
                nodes.append(node['coords'])
                node_indices.append(len(nodes))
            else:
                #added node
                node_indices.append(i)
        nurbs_indices.append(node_indices)
        
    cnt = 1
    fid.write('Begin Nodes\n')
    for node in nodes:
        fid.write(str(cnt))
        for i in range(0, len(node)):
            fid.write(' ')
            fid.write(node[i])
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End Nodes\n\n')
    
    fid.write('Begin Elements KinematicLinearGeo2dNURBS\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        mat = int(model_data['nurbs_data'][i]['material'])
        prop = model_data['materials_map'][mat]
        fid.write(str(prop))
        for node in nurbs_indices[i]:
            fid.write(' ')
            fid.write(str(node))
        fid.write('\n')
        cnt = cnt + 1
        model_data['nurbs_data'][i]['Node Indices'] = nurbs_indices[i]
    fid.write('End Elements\n\n')
    
    fid.write('Begin ElementalData NURBS_WEIGHTS\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        if 'Weights' in model_data['nurbs_data'][i]:
            size = len(model_data['nurbs_data'][i]['Weights'])
            fid.write('[' + str(size) + '] (')
            for j in range(0, size - 1):
                fid.write(model_data['nurbs_data'][i]['Weights'][j])
                fid.write(', ')
            fid.write(model_data['nurbs_data'][i]['Weights'][size - 1])
            fid.write(')\n')
            cnt = cnt + 1
        else:
            size = len(model_data['nurbs_data'][i]['Points'])
            fid.write('[' + str(size) + '] (')
            for j in range(0, size - 1):
                fid.write(str(1.0))
                fid.write(', ')
            fid.write(str(1.0))
            fid.write(')\n')
            cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.write('Begin ElementalData NURBS_KNOTS_1\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write('[' + str(len(model_data['nurbs_data'][i]['U Knots'])) + '] (')
        size = len(model_data['nurbs_data'][i]['U Knots'])
        for j in range(0, size - 1):
            fid.write(model_data['nurbs_data'][i]['U Knots'][j])
            fid.write(', ')
        fid.write(model_data['nurbs_data'][i]['U Knots'][size - 1])
        fid.write(')\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.write('Begin ElementalData NURBS_KNOTS_2\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write('[' + str(len(model_data['nurbs_data'][i]['V Knots'])) + '] (')
        size = len(model_data['nurbs_data'][i]['V Knots'])
        for j in range(0, size - 1):
            fid.write(model_data['nurbs_data'][i]['V Knots'][j])
            fid.write(', ')
        fid.write(model_data['nurbs_data'][i]['V Knots'][size - 1])
        fid.write(')\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.write('Begin ElementalData NURBS_DIMENSION_1\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write(model_data['nurbs_data'][i]['U Dim'])
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.write('Begin ElementalData NURBS_DIMENSION_2\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write(model_data['nurbs_data'][i]['V Dim'])
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.write('Begin ElementalData NURBS_DEGREE_1\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write(model_data['nurbs_data'][i]['U Degree'])
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.write('Begin ElementalData NURBS_DEGREE_2\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write(model_data['nurbs_data'][i]['V Degree'])
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    # TODO: add code to handle the number of division
    fid.write('Begin ElementalData NUM_DIVISION_1\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write(str(10))
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    # TODO: add code to handle the number of division
    fid.write('Begin ElementalData NUM_DIVISION_2\n')
    cnt = 1
    for i in range(0, len(model_data['nurbs_data'])):
        fid.write(str(cnt))
        fid.write(' ')
        fid.write(str(10))
        fid.write('\n')
        cnt = cnt + 1
    fid.write('End ElementalData\n\n')
    
    fid.close()

model_data = ReadRawData('rEpLaCeMeNtStRiNg_raw.dat')
GenerateMdpa('rEpLaCeMeNtStRiNg.mdpa', model_data)
print "Generate mdpa completed"

#print model_data

