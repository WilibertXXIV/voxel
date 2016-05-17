"""
    This script is designed to export a mass amount of MagicaVoxel .vox files
    to .obj. Unlike Magica's internal exporter, this exporter preserves the
    voxel vertices for easy manipulating in a 3d modeling program like Blender.
    
    This script also has an Obj optimizer (ObjViewer). The optimizer will
    attempt to merge adjacent same UV quads in a way that reduces vertices. This
    feature was added to aid in manually unwrapping the model.
    
    Notes:
        * ObjViewer is a misnomer, it can actually modify the file! Whoops...
            * Great pass! Sorry! What a save! What a save!
        * Both VoxelStruct and ObjViewer have .obj exporter code. This should
            probably be refactored to reduce code duplication.
        * There are a lot of floating point equality comparisons. They seem to
            work but it scares me a little.
"""
from operator import itemgetter
from functools import cmp_to_key

class ObjViewer:
    """ For reading OBJ files composed of axis aligned faces """
    def __init__(self):
        self.vertices = []
        self.faces = []
    def read(self, stream):
        """ Discards current model, loads a new one """
        self.vertices = []
        self.faces = []
        uvs = []
        normals = []
        for line in stream:
            # make sure there's no new line or trailing spaces
            l = line.strip().split(' ')
            lineType = l[0].strip()
            data = l[1:]
            if lineType == 'v':
                # vertex
                v = tuple(map(float, data))
                self.vertices.append(v)
            elif lineType == 'vt':
                # uv
                uvs.append( tuple(map(float, data)) )
            elif lineType == 'vn':
                # normal
                normals.append( tuple(map(float, data)) )
            elif lineType == 'f':
                # face (assume all verts/uvs/normals have been processed)
                faceVerts = []
                faceUvs = []
                faceNormals = []
                for v in data:
                    result = v.split('/')
                    # recall that everything is 1 indexed... #!@%!
                    faceVerts.append(self.vertices[int(result[0]) - 1])
                    if len(result) == 1:
                        continue
                    if result[1] != '':
                        # uvs may not be present, ex: 'f vert//normal ...'
                        faceUvs.append(uvs[int(result[1]) - 1])
                    if len(result) <= 2:
                        # don't continue if only vert and uv are present
                        continue
                    faceNormals.append(normals[int(result[2]) - 1])
                self.faces.append( ObjFace(faceVerts, faceUvs, faceNormals) )
    def optimize(self):
        """ Combine adjacent quads into bigger quads (finds a local max) """
        self.genNormals(False)
        # the edges combined with self.faces form an undirected graph
        adjacencyGraphEdges = self._buildAdjacencyGraph()
        groups = self._findGraphComponents(adjacencyGraphEdges)
        newFaces = []
        for group in groups:
            newFaces.extend(self._optimizeComponent(group, adjacencyGraphEdges))
        self._regen(newFaces)
    def _optimizeComponent(self, comp, edges):
        # if all faces are axis aligned, we can just splice out a coordinate and
        # use one simple algorithm to optimize the faces
        # recall that all faces in a component have the same UVs
        # (see facesAreAdjacent)
        normal = comp[0].normal
        if abs(normal[0]) == 1:
            # all verts have same x coordinate, so pass on (y, z)
            faces = self._optimizeSurface(comp, edges, 1, 2)
        elif abs(normal[1]) == 1:
            # all verts have same y coordinate, so pass on (x, z)
            faces = self._optimizeSurface(comp, edges, 0, 2)
        elif abs(normal[2]) == 1:
            # all verts have the same z coordinate, so pass on (x, y)
            faces = self._optimizeSurface(comp, edges, 0, 1)
        return faces
    def _regen(self, faces):
        self.vertices = set()
        self.faces = faces
        for f in self.faces:
            self.vertices.update(f.vertices)
        self.vertices = list(self.vertices)
    def _optimizeSurface(self, faces, edges, xPos, yPos):
        """ Optimize a 2d surface by breaking it into rectangles.
            See optimizeComponent for splicing explanation. xPos and yPos
            represent an index used to select coordintes from the vertex
            positions, which are represented by tuples.
        """
        pos = itemgetter(xPos, yPos)
        out = []
        rectRoot = None
        for face in faces:
            # TODO!
            pass
        return faces
    def _followStrip(self, face, edges, mask, output):
        adjacent = set(e for e in edges if e[0] is face)
        for f in adjacent:
            # TODO!
            pass
    def _buildAdjacencyGraph(self):
        """ Get the list of edges representing adjacent faces. """
        # a list of edges between adjacent face tuple(face_a, face_b)
        self.edges = []
        # build the list of edges in the graph
        for root in self.faces:
            for face in self.faces:
                if face is root:
                    continue
                if self.facesAreAdjacent(root, face):
                    # the other edge will happen somewhere else in the iteration
                    # (i.e., the relation isAdjacent is symmetric)
                    self.edges.append((root, face))
        return self.edges
    def _findGraphComponents(self, edges):
        """ Get the list of connected components from a list of graph edges.
            The list will contain lists containing the edges of the graph.
            
            The result is stored in groups.
        """
        groups = []
        visited = dict((f, False) for f in self.faces)
        for face in self.faces:
            # if the face hasn't been visited, it is not in any found components
            if not visited[face]:
                g = []
                self._visitGraphNodes(face, edges, visited, g)
                # there is only a new component if face has not been visited yet
                groups.append(g)
        return groups
    def _visitGraphNodes(self, node, edges, visited, component):
        # visit every component connected to this one
        for edge in edges:
            # for all x in nodes, (node, x) and (x, node) should be in edges!
            # therefore we don't have to check for "edge[1] is node"
            if edge[0] is node and not visited[edge[1]]:
                assert edge[1] is not node
                # mark the other node as visited
                visited[edge[1]] = True
                component.append(edge[1])
                # visit all of that nodes connected nodes
                self._visitGraphNodes(edge[1], edges, visited, component)
    def facesAreAdjacent(self, a, b):
        """ Adjacent is defined as save normal, uv, and a shared edge.
            This isn't entirely intuitive (i.e., corner faces are not adjacent)
            but this definition fits the problem domain.
        """
        # note: None is == None, this shouldn't matter
        if a.uv != b.uv:
            return False
        if a.normal != b.normal:
            return False
        # to be adjacent, two faces must share an edge
        # use == and not identity in case the edge split was used
        shared = 0
        for vert_a in a.vertices:
            for vert_b in b.vertices:
                if vert_a == vert_b:
                    shared += 1
                # hooray we have found a shared edge (or a degenerate case...)
                if shared == 2:
                    return True
        return False
    def genNormals(self, overwrite=False):
        # compute CCW normal if it doesn't exist
        for face in self.faces:
            if overwrite or face.normal is None:
                side_a = (face.vertices[1][0] - face.vertices[0][0],
                          face.vertices[1][1] - face.vertices[0][1],
                          face.vertices[1][2] - face.vertices[0][2])
                side_b = (face.vertices[-1][0] - face.vertices[0][0],
                          face.vertices[-1][1] - face.vertices[0][1],
                          face.vertices[-1][2] - face.vertices[0][2])
                # compute the cross product
                face.normal = (side_a[1]*side_b[2] - side_a[2]*side_b[1],
                               side_a[2]*side_b[0] - side_a[0]*side_b[2],
                               side_a[0]*side_b[1] - side_a[1]*side_b[0])
    def exportObj(self, stream):
        # gather some of the needed information
        normals = set()
        uvs = set()
        for f in self.faces:
            if f.uv is not None:
                uvs.add(f.uv)
            if f.normal is not None:
                normals.add(f.normal)
        # convert these to lists because we need to get their index later
        normals = list(normals)
        uvs = list(uvs)
        # write to the file
        stream.write('# shivshank\'s .obj optimizer\n')
        stream.write('\n')
        if len(normals) > 0:
            stream.write('# normals\n')
            for n in normals:
                stream.write('vn ' + ' '.join(list(map(str, n))) + '\n')
            stream.write('\n')
        if len(uvs) > 0:
            stream.write('# texcoords\n')
            for i in uvs:
                stream.write('vt ' + ' '.join(list(map(str, i))) + '\n')
            stream.write('\n')
        stream.write('# verts\n')
        for v in self.vertices:
            stream.write('v ' + ' '.join(list(map(str, v))) + '\n')
        stream.write('\n')
        stream.write('# faces\n')
        for f in self.faces:
            # recall that OBJ files are 1 indexed
            if f.normal is not None:
                n = 1 + normals.index(f.normal)
            else:
                n = ''
            if f.uv is not None:
                uv = 1 + normals.index(f.uv)
            else:
                uv = ''
            # this used to be a one liner ;)
            fLine = []
            for vert in f.vertices:
                # I think this is going to be a speed bottleneck!
                v = 1 + self.vertices.index(vert)
                fLine.append(str(v) + '/' + str(uv) + '/' + str(n))
            fLine = ' '.join(fLine)
            stream.write('f ' + fLine + '\n')
        stream.write('\n')
        stream.write('\n')

class ObjFace:
    """ An axis aligned quad. """
    def __init__(self, verts, uvs, normals):
        self.vertices = verts
        assert len(verts) == 4, "only quads are supported"
        for i in normals:
            assert normals[0] == i, "face must be axis aligned (orthogonal normals)"
        self.normal = normals[0] if len(normals) > 0 else None
        for i in uvs:
            assert uvs[0] == i, "face must be axis aligned (orthogonal normals)"
        self.uv = uvs[0] if len(uvs) > 0 else None
        self.center = (
            sum(i[0] for i in self.vertices)/4,
            sum(i[1] for i in self.vertices)/4,
            sum(i[2] for i in self.vertices)/4
        )

class VoxelStruct:
    """ Describes a voxel object
    """
    def __init__(self, voxelList):
        # accessed by (x, y, z) {hopefully python can efficiently hash this}
        self.voxels = {}
        # the bounds are defined by a minimum point and a maximum point
        self.origin = (0, 0, 0)
        self.maximum = (1, 1, 1)
        self.colorIndices = set()
        self._process(voxelList)
    def _process(self, voxels):
        for voxel in voxels:
            self.voxels[(voxel.x, voxel.y, voxel.z)] = voxel
            self.colorIndices.add(voxel.colorIndex)
            # update the bounds
            self.origin = (
                min(self.origin[0], voxel.x),
                min(self.origin[1], voxel.y),
                min(self.origin[2], voxel.z)
            )
            self.maximum = (
                max(self.maximum[0], voxel.x),
                max(self.maximum[1], voxel.y),
                max(self.maximum[2], voxel.z)
            )
    def exportObj(self, stream):
        # produces an obj with normals, UVs, and vertex coordinates
        # note: texcoords are based on MagicaVoxel's texturing scheme!
        #   meaning a color index of 0 translates to pixel[255]
        #   and color index [1:256] -> pixel[0:255]
        uv = []
        for i in self.colorIndices:
            if i == 0:
                continue
            # all UV v coordinates are 0.5, along middle of the texture
            uv.append(((i - 1) / 256 + 1/512, 0.5))
        # populate the normals
        normals = [
            (-1, 0, 0), # left = 0
            (1, 0, 0),  # right = 1
            (0, 0, 1),  # top = 2
            (0, 0, -1), # bottom = 3
            (0, -1, 0), # front = 4
            (0, 1, 0)   # back = 5
        ]
        verts = []
        faces = []
        for pos, voxel in self.voxels.items():
            vertIndex = len(verts)
            newVerts = list(self._getObjVerts(voxel))
            verts.extend(newVerts)
            faces.extend(self._getObjFaces(voxel, vertIndex, newVerts, uv))
        # write the result!
        stream.write('# shivshank\'s .vox exporter\n')
        stream.write('\n')
        stream.write('# normals\n')
        for i in normals:
            stream.write('vn ' + ' '.join(list(map(str, i))) + '\n')
        stream.write('\n')
        stream.write('# texcoords\n')
        for i in uv:
            stream.write('vt ' + ' '.join(list(map(str, i))) + '\n')
        stream.write('\n')
        stream.write('# verts\n')
        for i in verts:
            stream.write('v ' + ' '.join(list(map(str, i))) + '\n')
        stream.write('\n')
        # n.b., faces actually holds strings to avoid another data structure
        stream.write('# faces\n')
        for i in faces:
            stream.write('f ' + i + '\n')
        stream.write('\n')
        stream.write('\n')
    def _getObjFaces(self, voxel, baseVertIndex, verts, uv):
        if voxel.colorIndex == 0:
            # do nothing if this is an empty voxel
            return
        # face and vertex access should really be combined as to reduce
        # redundant computations
        exposed = self._objExposed(voxel)
        faces = []
        for side in exposed:
            if side == 0:
                f = [
                    1 + verts.index((voxel.x, voxel.y + 1, voxel.z + 1)),
                    1 + verts.index((voxel.x, voxel.y + 1, voxel.z)),
                    1 + verts.index((voxel.x, voxel.y, voxel.z)),
                    1 + verts.index((voxel.x, voxel.y, voxel.z + 1))
                ]
            elif side == 1:
                f = [
                    1 + verts.index((voxel.x + 1, voxel.y, voxel.z + 1)),
                    1 + verts.index((voxel.x + 1, voxel.y, voxel.z)),
                    1 + verts.index((voxel.x + 1, voxel.y + 1, voxel.z)),
                    1 + verts.index((voxel.x + 1, voxel.y + 1, voxel.z + 1))
                ]
            elif side == 2:
                f = [
                    1 + verts.index((voxel.x, voxel.y + 1, voxel.z + 1)),
                    1 + verts.index((voxel.x, voxel.y, voxel.z + 1)),
                    1 + verts.index((voxel.x + 1, voxel.y, voxel.z + 1)),
                    1 + verts.index((voxel.x + 1, voxel.y + 1, voxel.z + 1))
                ]
            elif side == 3:
                f = [
                    1 + verts.index((voxel.x, voxel.y, voxel.z)),
                    1 + verts.index((voxel.x, voxel.y + 1, voxel.z)),
                    1 + verts.index((voxel.x + 1, voxel.y + 1, voxel.z)),
                    1 + verts.index((voxel.x + 1, voxel.y, voxel.z))
                ]
            elif side == 4:
                f = [
                    1 + verts.index((voxel.x, voxel.y, voxel.z + 1)),
                    1 + verts.index((voxel.x, voxel.y, voxel.z)),
                    1 + verts.index((voxel.x + 1, voxel.y, voxel.z)),
                    1 + verts.index((voxel.x + 1, voxel.y, voxel.z + 1))
                ]
            elif side == 5:
                f = [
                    1 + verts.index((voxel.x + 1, voxel.y + 1, voxel.z + 1)),
                    1 + verts.index((voxel.x + 1, voxel.y + 1, voxel.z)),
                    1 + verts.index((voxel.x, voxel.y + 1, voxel.z)),
                    1 + verts.index((voxel.x, voxel.y + 1, voxel.z + 1))
                ]
            else:
                # this should never happen!
                assert False
            # the normal is at index 0, so in obj file it's index 1
            n = str(1 + side)
             # floating point error could make this tricky -.-
            u = uv.index( ((voxel.colorIndex - 1)/256 + 1/512, 0.5) )
            u = str(1 + u)
            faces.append(
                str(f[0] + baseVertIndex) + '/' + u + '/' + n + ' ' +
                str(f[1] + baseVertIndex) + '/' + u + '/' + n + ' ' +
                str(f[2] + baseVertIndex) + '/' + u + '/' + n + ' ' +
                str(f[3] + baseVertIndex) + '/' + u + '/' + n
            )
        return faces
    def _getObjVerts(self, voxel):
        # MagicaVoxel does -.5 to +.5 for each cube, we'll do 0.0 to 1.0 ;)
        exposed = self._objExposed(voxel)
        verts = set()
        if 0 in exposed:
            verts.add((voxel.x, voxel.y, voxel.z))
            verts.add((voxel.x, voxel.y, voxel.z + 1))
            verts.add((voxel.x, voxel.y + 1, voxel.z))
            verts.add((voxel.x, voxel.y + 1, voxel.z + 1))
        if 1 in exposed:
            verts.add((voxel.x + 1, voxel.y, voxel.z))
            verts.add((voxel.x + 1, voxel.y, voxel.z + 1))
            verts.add((voxel.x + 1, voxel.y + 1, voxel.z))
            verts.add((voxel.x + 1, voxel.y + 1, voxel.z + 1))
        if 2 in exposed:
            verts.add((voxel.x, voxel.y, voxel.z + 1))
            verts.add((voxel.x, voxel.y + 1, voxel.z + 1))
            verts.add((voxel.x + 1, voxel.y, voxel.z + 1))
            verts.add((voxel.x + 1, voxel.y + 1, voxel.z + 1))
        if 3 in exposed:
            verts.add((voxel.x, voxel.y, voxel.z))
            verts.add((voxel.x, voxel.y + 1, voxel.z))
            verts.add((voxel.x + 1, voxel.y, voxel.z))
            verts.add((voxel.x + 1, voxel.y + 1, voxel.z))
        if 4 in exposed:
            verts.add((voxel.x, voxel.y, voxel.z + 1))
            verts.add((voxel.x, voxel.y, voxel.z))
            verts.add((voxel.x + 1, voxel.y, voxel.z))
            verts.add((voxel.x + 1, voxel.y, voxel.z + 1))
        if 5 in exposed:
            verts.add((voxel.x + 1, voxel.y + 1, voxel.z + 1))
            verts.add((voxel.x + 1, voxel.y + 1, voxel.z))
            verts.add((voxel.x, voxel.y + 1, voxel.z))
            verts.add((voxel.x, voxel.y + 1, voxel.z + 1))
        return verts
    def _objExposed(self, voxel):
        """ get the sick truth about these voxels' dirty secrets
        """
        res = set(i for i in range(6))
        # check left 0
        side = self.voxels.get((voxel.x - 1, voxel.y, voxel.z))
        if side is not None and side.colorIndex != 0:
            res.remove(0)
        # check right 1
        side = self.voxels.get((voxel.x + 1, voxel.y, voxel.z))
        if side is not None and side.colorIndex != 0:
            res.remove(1)
        # check top 2
        side = self.voxels.get((voxel.x, voxel.y, voxel.z + 1))
        if side is not None and side.colorIndex != 0:
            res.remove(2)
        # check bottom 3
        side = self.voxels.get((voxel.x, voxel.y, voxel.z - 1))
        if side is not None and side.colorIndex != 0:
            res.remove(3)
        # check front 4
        side = self.voxels.get((voxel.x, voxel.y - 1, voxel.z))
        if side is not None and side.colorIndex != 0:
            res.remove(4)
        # check back 5
        side = self.voxels.get((voxel.x, voxel.y + 1, voxel.z))
        if side is not None and side.colorIndex != 0:
            res.remove(5)
        return res

class Voxel:
    def __init__(self, x, y, z, colorIndex):
        self.x = x
        self.y = y
        self.z = z
        self.colorIndex = colorIndex

def decodeVox(path):
    # in theory this could elegantly be many functions and classes
    # but this is such a simple file format...
    voxels = []
    with open(path, mode='rb') as file:
        magic = file.read(4)
        if magic != b'VOX ':
            print('magic number is', magic)
            if userAborts('This does not appear to be a VOX file. Abort?'):
                return
        # the file appears to use little endian consistent with RIFF
        version = int.from_bytes(file.read(4), byteorder='little')
        if version != 150:
            userAborts('Only version 150 is supported; this file: '
                      + str(version) + '. Abort?')
        print('\treading main chunk')
        mainHeader = readChunkHeader(file)
        if mainHeader['id'] != b'MAIN':
            print('chunk id:', mainId)
            if userAborts('Did not find the main chunk. Abort?'):
                return
        # main should be empty
        assert mainHeader['size'] == 0
        nextHeader = readChunkHeader(file)
        while nextHeader['id'] != b'XYZI':
            # we don't need anything from the size or palette header!
            # we'll figure out (minimum) bounds later from the voxel data
            # and we know what palette we're using
            nextHeader = readChunkHeader(file)
        voxelHeader = nextHeader
        assert voxelHeader['id'] == b'XYZI'
        assert voxelHeader['childrenSize'] == 0
        seekPos = file.tell()
        totalVoxels = int.from_bytes(file.read(4), byteorder='little')
        ### READ THE VOXELS ###
        for i in range(totalVoxels):
            # n.b., byte order should be irrelevant since these are all 1 byte
            x = int.from_bytes(file.read(1), byteorder='little')
            y = int.from_bytes(file.read(1), byteorder='little')
            z = int.from_bytes(file.read(1), byteorder='little')
            color = int.from_bytes(file.read(1), byteorder='little')
            voxels.append(Voxel(x, y, z, color))
        # (there may be more chunks after this but we don't need them!)
        # assert that we've read the entire voxel chunk
        assert file.tell() - seekPos == voxelHeader['size']
        print('\tdone reading voxel data;', totalVoxels , 'voxels read ;D')
    return voxels

def readChunkHeader(buffer):
    id = buffer.read(4)
    size = int.from_bytes(buffer.read(4), byteorder='little')
    childrenSize = int.from_bytes(buffer.read(4), byteorder='little')
    return {
        'id': id, 'size': size, 'childrenSize': childrenSize
    }

def userAborts(msg):
    print(msg + ' (y/n)')
    u = input()
    if u.startswith('n'):
        return False
    
    return True

if __name__ == '__main__':
    import os, os.path
    from glob import glob
    
    print('Enter an output path:')
    u = input('> ').strip()
    while not os.path.exists(u):
        print('That path does not exist.')
        print('Enter an output path:')
        u = input('> ').strip()
    outPath = os.path.abspath(u)
    
    print('Are we optimizing? (y/n)')
    u = input('> ').strip()
    if u.startswith('y'):
        optimizing = True
    else:
        optimizing = False

    try:
        while True:
            print('Enter glob of export files (\'exit\' or blank to quit):')
            u = input('> ').strip()
            if u == 'exit' or u == '':
                break
            u = glob(u)
            for f in u:
                print('reading VOX file', f)
                vox = decodeVox(f)
                if vox is None:
                    continue
                res = VoxelStruct(vox)
                out = os.path.join(outPath,
                              os.path.splitext(os.path.basename(f))[0] + '.obj')
                print('exporting VOX to OBJ at path', out)
                with open(out, mode='w') as file:
                    res.exportObj(file)
                if optimizing:
                    print('optimizing OBJ at path', out)
                    with open(out, mode='r') as file:
                        opti = ObjViewer()
                        opti.read(file)
                    opti.optimize()
                    print('exporting optimized OBJ to', out)
                    with open(out, mode='w') as file:
                        opti.exportObj(file)
                    
    except KeyboardInterrupt:
        pass