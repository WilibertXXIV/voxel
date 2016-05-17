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
                v = tuple(map(int, data))
                self.vertices.append(v)
            elif lineType == 'vt':
                # uv
                pass
            elif lineType == 'vn':
                # normal
                pass
            elif lineType == 'f':
                # face (assume all verts/uvs/normals have been processed)
                faceVerts = []
                faceUvs = []
                faceNormals = []
                for v in data:
                    result = v.split('/')
                    faceVerts.append(self.vertices[result[0]])
                    faceUvs.append(uvs[result[1]])
                    faceNormals.append(normals[self.vertices[result[2]]])
                self.faces.append( ObjFace(faceVerts, faceUvs, faceNormals) )
    def optimize(self):
        """ Combine adjacent quads into bigger quads (finds a local max) """
        pass

class ObjFace:
    def __init__(self, verts):
        pass

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
    except KeyboardInterrupt:
        pass