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
        * TODO: ObjViewer currently is not used and does nothing useful. It is
            currently being refactored.
        * There are a few floating point equality comparisons. They seem to work
            but it scares me a little.
        * TODO: use constants instead of magic numbers (as defined in AAQuad),
                (i.e., ..., 2 -> AAQuad.TOP, ...)
        * A lot of assertions should probably be exceptions since they are
            error checking user input (this sounds really bad now that I've put
            it on paper...). So don't run in optimized mode (who does that
            anyways?).
"""
from operator import itemgetter
from functools import cmp_to_key

class ObjViewer:
    """ For reading OBJ files composed of axis aligned faces """
    def __init__(self):
        self.vertices = []
        self.faces = []
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

class AAQuad:
    """ A solid colored axis aligned quad. """
    normals = [
        (-1, 0, 0), # left = 0
        (1, 0, 0),  # right = 1
        (0, 0, 1),  # top = 2
        (0, 0, -1), # bottom = 3
        (0, -1, 0), # front = 4
        (0, 1, 0)   # back = 5
    ]
    LEFT = 0
    RIGHT = 1
    TOP = 2
    BOTTOM = 3
    FRONT = 4
    BACK = 5
    def __init__(self, verts, uv=None, normal=None):
        assert len(verts) == 4, "face must be a quad"
        self.vertices = verts
        self.uv = uv
        self.normal = normal
    def __str__(self):
        s = []
        for i in self.vertices:
            s.append( str(i) + '/' + str(self.uv) + '/' + str(self.normal))
        return 'f ' + ' '.join(s)

class ObjFace:
    """ An arbitrary geometry face """
    def __init__(self, verts, uvs=None, normals=None):
        self.vertices = verts
        assert len(verts) in (3, 4), "only quads and tris are supported"
        self.normals = normals
        self.uvs = uvs
    def toAAQuad(self):
        q = AAQuad(self.vertices)
        if self.normals is not None and len(self.normals) > 0:
            for i in self.normals:
                assert self.normals[0] == i, \
                    "face must be axis aligned (orthogonal normals)"
            q.normal = self.normals[0]
        if self.uvs is not None and len(self.uvs) > 0:
            for i in self.uvs:
                assert self.uvs[0] == i, \
                    "face must be axis aligned (orthogonal)"
            q.uv = self.uvs[0]
        return q

class VoxelStruct:
    """ Describes a voxel object
    """
    def __init__(self, voxelList):
        # a dict is probably the best way to go about this
        # (as a trade off between performance and code complexity)
        # accessed by (x, y, z) {hopefully python can efficiently hash this}
        self.voxels = {}
        # the bounds are defined by a minimum point and a maximum point
        self.origin = (0, 0, 0)
        self.maximum = (0, 0, 0)
        self.colorIndices = set()
        self._process(voxelList)
    def _process(self, voxels):
        self.voxels = {}
        self.origin = (0, 0, 0)
        self.maximum = (0, 0, 0)
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
    def zeroOrigin(self):
        """ Translate the model so that it's origin is at 0, 0, 0 """
        result = {}
        xOff, yOff, zOff = self.origin
        for key, voxel in self.voxels.items():
            result[(voxel.x-xOff, voxel.y-yOff, voxel.z-zOff)] = \
                Voxel(voxel.x-xOff, voxel.y-yOff, voxel.z-zOff,
                      voxel.colorIndex)
        self.voxels = result
        self.origin = (0, 0, 0)
        self.maximum = (self.maximum[0] - xOff,
                        self.maximum[1] - yOff,
                        self.maximum[2] - zOff)
    def toQuads(self):
        """ --> a list of AAQuads """
        faces = []
        for pos, voxel in self.voxels.items():
            faces.extend(self._getObjFaces(voxel))
        return faces
    def _getObjFaces(self, voxel):
        if voxel.colorIndex == 0:
            # do nothing if this is an empty voxel
            # n.b., I do not know if this ever can happen.
            return []
        exposed = self._objExposed(voxel)
        faces = []
        for side in exposed:
            if side == 0:
                f = self._getLeftSide(voxel)
            elif side == 1:
                f = self._getRightSide(voxel)
            elif side == 2:
                f = self._getTopSide(voxel)
            elif side == 3:
                f = self._getBottomSide(voxel)
            elif side == 4:
                f = self._getFrontSide(voxel)
            elif side == 5:
                f = self._getBackSide(voxel)
            else:
                assert False, 'a weird number ended up in the exposed mask/set'
            n = AAQuad.normals[side]
            # note: texcoords are based on MagicaVoxel's texturing scheme!
            #   meaning a color index of 0 translates to pixel[255]
            #   and color index [1:256] -> pixel[0:255]
            u = ((voxel.colorIndex - 1)/256 + 1/512, 0.5)
            faces.append(
                # this is most definitely not "fun"
                AAQuad(f, u, n)
            )
        return faces
    # MagicaVoxel does -.5 to +.5 for each cube, we'll do 0.0 to 1.0 ;)
    def _getLeftSide(self, voxel):
        return [
            (voxel.x, voxel.y + 1, voxel.z + 1),
            (voxel.x, voxel.y + 1, voxel.z),
            (voxel.x, voxel.y, voxel.z),
            (voxel.x, voxel.y, voxel.z + 1)
        ]
    def _getRightSide(self, voxel):
        return [
            (voxel.x + 1, voxel.y, voxel.z + 1),
            (voxel.x + 1, voxel.y, voxel.z),
            (voxel.x + 1, voxel.y + 1, voxel.z),
            (voxel.x + 1, voxel.y + 1, voxel.z + 1)
        ]
    def _getTopSide(self, voxel):
        return [
            (voxel.x, voxel.y + 1, voxel.z + 1),
            (voxel.x, voxel.y, voxel.z + 1),
            (voxel.x + 1, voxel.y, voxel.z + 1),
            (voxel.x + 1, voxel.y + 1, voxel.z + 1)
        ]
    def _getBottomSide(self, voxel):
        return [
            (voxel.x, voxel.y, voxel.z),
            (voxel.x, voxel.y + 1, voxel.z),
            (voxel.x + 1, voxel.y + 1, voxel.z),
            (voxel.x + 1, voxel.y, voxel.z)
        ]
    def _getFrontSide(self, voxel):
        return [
            (voxel.x, voxel.y, voxel.z + 1),
            (voxel.x, voxel.y, voxel.z),
            (voxel.x + 1, voxel.y, voxel.z),
            (voxel.x + 1, voxel.y, voxel.z + 1)
        ]
    def _getBackSide(self, voxel):
        return [
            (voxel.x + 1, voxel.y + 1, voxel.z + 1),
            (voxel.x + 1, voxel.y + 1, voxel.z),
            (voxel.x, voxel.y + 1, voxel.z),
            (voxel.x, voxel.y + 1, voxel.z + 1)
        ]
    def _objExposed(self, voxel):
        """ --> a set of [0, 6) representing which voxel faces are shown
            for the meaning of 0-5, see AAQuad.normals
            get the sick truth about these voxels' dirty secrets...
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

def importObj(stream):
    vertices = []
    faces = []
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
            vertices.append(v)
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
                print(result)
                # recall that everything is 1 indexed...
                faceVerts.append(vertices[int(result[0]) - 1])
                if len(result) == 1:
                    continue # there is only a vertex index
                if result[1] != '':
                    # uvs may not be present, ex: 'f vert//normal ...'
                    faceUvs.append(uvs[int(result[1]) - 1])
                if len(result) <= 2:
                    # don't continue if only vert and uv are present
                    continue
                faceNormals.append(normals[int(result[2]) - 1])
            faces.append( ObjFace(faceVerts, faceUvs, faceNormals) )
        else:
            # there could be material specs, smoothing, or comments... ignore!
            pass
    return faces

def exportObj(stream, aaQuads):
    # gather some of the needed information
    faces = aaQuads
    normals = set()
    uvs = set()
    for f in faces:
        if f.uv is not None:
            uvs.add(f.uv)
        if f.normal is not None:
            normals.add(f.normal)
    # convert these to lists because we need to get their index later
    normals = list(normals)
    uvs = list(uvs)
    # we will build a list of vertices as we go and then write everything
    # in bulk, disadvantage that MANY verts will be duplicated in the OBJ file
    fLines = []
    vertices = []
    indexOffset = 0
    for f in faces:
        # recall that OBJ files are 1 indexed
        n = 1 + normals.index(f.normal) if f.normal is not None else ''
        uv = 1 + normals.index(f.uv) if f.uv is not None else ''
        # this used to be a one liner ;)
        fLine = ['f']
        for i, vert in enumerate(f.vertices):
            # for each vertex of this face
            v = 1 + indexOffset + f.vertices.index(vert)
            fLine.append(str(v) + '/' + str(uv) + '/' + str(n))
        vertices.extend(f.vertices)
        indexOffset += len(f.vertices)
        fLines.append(' '.join(fLine) + '\n')
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
    # output the vertices and faces
    stream.write('# verts\n')
    for v in vertices:
        stream.write('v ' + ' '.join(list(map(str, v))) + '\n')
    stream.write('\n')
    stream.write('# faces\n')
    for i in fLines:
        stream.write(i)
    stream.write('\n')
    stream.write('\n')

def importVox(file):
    """ --> a VoxelStruct from this .vox file stream """
    # in theory this could elegantly be many functions and classes
    # but this is such a simple file format...
    # refactor: ? should probably find a better exception type than value error
    voxels = []
    magic = file.read(4)
    if magic != b'VOX ':
        print('magic number is', magic)
        if userAborts('This does not appear to be a VOX file. Abort?'):
            raise ValueError("Invalid magic number")
    # the file appears to use little endian consistent with RIFF
    version = int.from_bytes(file.read(4), byteorder='little')
    if version != 150:
        if userAborts('Only version 150 is supported; this file: '
                      + str(version) + '. Abort?'):
            raise ValueError("Invalid file version")
    mainHeader = readChunkHeader(file)
    if mainHeader['id'] != b'MAIN':
        print('chunk id:', mainId)
        if userAborts('Did not find the main chunk. Abort?'):
            raise ValueError("Did not find main VOX chunk. ")
    assert mainHeader['size'] == 0, "main chunk should have size 0"
    # we don't need anything from the size or palette header!
    # : we can figure out (minimum) bounds later from the voxel data
    # : we only need UVs from voxel data; user can export palette elsewhere
    nextHeader = readChunkHeader(file)
    while nextHeader['id'] != b'XYZI':
        # skip the contents of this header and its children, read the next one
        file.read(nextHeader['size'] + nextHeader['childrenSize'])
        nextHeader = readChunkHeader(file)
    voxelHeader = nextHeader
    assert voxelHeader['id'] == b'XYZI', 'this should be literally impossible'
    assert voxelHeader['childrenSize'] == 0, 'why voxel chunk have children?'
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
    # assert that we've read the entire voxel chunk
    assert file.tell() - seekPos == voxelHeader['size']
    # (there may be more chunks after this but we don't need them!)
    #print('\tdone reading voxel data;', totalVoxels , 'voxels read ;D')
    vs = VoxelStruct(voxels)
    return vs
    
def readChunkHeader(buffer):
    id = buffer.read(4)
    if id == b'':
        raise ValueError("Unexpected EOF, expected chunk header")
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

def exportAll():
    """ Uses a file to automatically export a bunch of files!
        See this function for details on the what the file looks like.
    """
    import os, os.path

    with open('exporter.txt', mode='r') as file:
        # use this as a file "spec"
        fromSource = os.path.abspath(file.readline().strip())
        toExportDir = os.path.abspath(file.readline().strip())
        optimizing = file.readline()
        if optimizing.lower() == 'true':
            optimizing = True
        else:
            optimizing = False

    print('exporting vox files under', fromSource)
    print('\tto directory', toExportDir)
    print('\toptimizing?', optimizing)
    print()
    
    # export EVERYTHING (.vox) walking the directory structure
    for p, dirList, fileList in os.walk(fromSource):
        pathDiff = os.path.relpath(p, start=fromSource)
        outDir = os.path.join(toExportDir, pathDiff)
        for fileName in fileList:
            if os.path.splitext(fileName)[1] != '.vox':
                # only take vox files
                print('\tignored', fileName)
                continue
            # read the voxel file
            with open(os.path.join(p, fileName), mode='r') as file:
                print('\texporting', fileName)
                vox = decodeVox(file)
            # mirror the directory structure in the export folder
            if not os.path.exists(outDir):
                os.makedirs(outDir)
                print('\tcreated directory', outDir)
            # always export a non-optimized version
            objName = os.path.splitext(fileName)[0]
            with open(os.path.join(outDir, objName + '.obj'), mode='w') as file:
                # convert it to AAQuads
                exportObj(file, vox.toQuads())
            if optimizing:
                # TODO: export the optimized version
                pass
                """obj = ObjViewer()
                with open(os.path.join(outDir, objName + '.obj'),
                          mode='w') as file:
                    obj.read(file)
                obj.optimize()
                with open(os.path.join(outDir, objName + '.opti.obj'),
                          mode='w') as file:
                    obj.exportObj(file)"""

def byPrompt():
    import os, os.path
    from glob import glob

    print('Enter an output path:')
    u = input('> ').strip()
    while not os.path.exists(u):
        print('That path does not exist.')
        print('Enter an output path:')
        u = input('> ').strip()
    outRoot = os.path.abspath(u)
    
    print('Are we optimizing? (y/n)')
    u = input('> ').strip()
    # this could be a one liner but I think it's easier to read this way
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
                with open(f, mode='r') as file:
                    try:
                        vox = decodeVox(file)
                    except ValueError:
                        print('\tfile reading aborted')
                        continue
                outFile = os.path.splitext(os.path.basename(f))[0]
                outPath = os.path.join(outRoot, outFile)
                print('exporting VOX to OBJ at path', outPath)
                with open(outPath, mode='w') as file:
                    exportObj(file, vox.toQuads())
                if optimizing:
                    # TODO
                    pass
                    """
                    print('optimizing OBJ at path', out)
                    with open(out, mode='r') as file:
                        opti = ObjViewer()
                        opti.read(file)
                    opti.optimize()
                    print('exporting optimized OBJ to', out)
                    with open(out, mode='w') as file:
                        opti.exportObj(file)"""
    except KeyboardInterrupt:
        pass

if __name__ == "__main__":
    try:
        exportAll()
    except OSError:
        print('No instruction file found, falling back to prompt.')
        byPrompt()