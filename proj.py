# Alim Dhanani, 101156584

GEOMPROC = "../GeomProc"
import sys
sys.path.append(GEOMPROC)
import geomproc
import numpy as np
from scipy.spatial import Delaunay
import time
AVERAGE_TYPES = ["voronoi", "mean", "corner"]
VORONOI_AVERAGE = 0
MEAN_AVERAGE = 1
CORNER_AVERAGE = 2

RED = np.array([1,0,0])
GREEN = np.array([0,1,0])
BLUE = np.array([0,0,1])

class Grid:
    def __init__(self, num_cubes, dimension=3, average_type=0):
        self.start = -1
        self.end = 1
        self.num_cubes = num_cubes
        self.dimension = dimension
        self.average_type = average_type
        self.cube_size = (self.end-self.start) / self.num_cubes

        self.next_vertex_id = VORONOI_AVERAGE
        self.vertices_in_cell = {}
        self.vertex_id = {}

    def getCell(self, point) -> np.array:
        point = np.asarray(point)
        if (len(point) != self.dimension):
            raise Exception("Invalid argument to getCell(): invalid point dimensionality")
        output = np.zeros(self.dimension)

        for i in range(len(point)):
            coord = point[i]
            output[i] = ((coord + 1) // self.cube_size)

        if (output == 0).all():
            raise Exception(f"no cell found for point {point}")
        return output
    
    def getCellCorner(self, cell) -> np.array:
        cell = np.asarray(cell)
        if (len(cell) != self.dimension):
            raise Exception("Invalid argument to getCellCoordinate(): invalid cell dimensionality")
        corner = np.zeros(self.dimension)

        for i in range(len(cell)):
            cube = cell[i]
            corner[i] = (cube * self.cube_size) - 1

        return corner
    

    def calculateMergedVertices(self, mesh):
        '''
        Calculating each cell's output vertex and vertex normal
        '''
        if (mesh.vnormal.shape[0] == 0):
            mesh.compute_vertex_and_face_normals()

        # initialize new-vertices as a zero-filled 2D array with dimensions ("next-vertex-id", mesh.dimension)
        new_vertices = np.zeros((self.next_vertex_id, self.dimension), dtype=float)

        # initialize new-vnormals as a zero-filled 2D array with dimensions ("next-vertex-id", mesh.dimension)
        new_vnormals = np.zeros((self.next_vertex_id, self.dimension), dtype=float)

        # get voronoi areas of all vertices in original mesh
        vareas = calculateMeshVoronoi(mesh)

        for cell in self.vertex_id:
            # for each cell c_i in vertex-id

            # retrieve the new vertex id
            new_id = self.vertex_id[cell]

            # initialize cell-vareas
            cell_vareas = vareas[self.vertices_in_cell[cell]]

            # total-varea = sum(cell-varea)
            total_varea = np.sum(cell_vareas)

            for j in range(len(self.vertices_in_cell[cell])):
            # for each vertex $v_j$ in vertices-in-cell[$hash(c_i)$]
                v_j = self.vertices_in_cell[cell][j]

                if (self.average_type == VORONOI_AVERAGE):
                    # new-vertices[new-id] = weighted mean of every v_j's coordinates, weighted by v_j's voronoi area
                    new_vertices[new_id, :] += (cell_vareas[j] / total_varea) * mesh.vertex[v_j]

                    new_vnormals[new_id, :] += (cell_vareas[j] / total_varea) * mesh.vnormal[v_j]

                elif (self.average_type == MEAN_AVERAGE):
                    new_vertices[new_id, :] += mesh.vertex[v_j] / len(self.vertices_in_cell[cell])
                    new_vnormals[new_id, :] += mesh.vnormal[v_j] / len(self.vertices_in_cell[cell])

            if (self.average_type == CORNER_AVERAGE):
                cell_pos = self.getCell(mesh.vertex[self.vertices_in_cell[cell][0]])
                new_vertices[new_id, :] = self.getCellCorner(cell_pos)

        return new_vertices, new_vnormals

    def simplifyMesh(self, mesh:geomproc.mesh):
        mesh.compute_vertex_and_face_normals()
        new_triangles = []
        face_colors = []
        new_triangles_dict = {}
        for t in range(len(mesh.face)):
            # for each triangle t_{in} in original mesh

            t_in = mesh.face[t]
            cells = np.zeros((3, self.dimension), dtype=np.int_)
            hashes = np.zeros(3)

            for i in range(3):
                # for each vertex index v_i in triangle t_{in}
    
                # find cell c = (c_1, c_2, c_3) using math of v_i's coordinates
                (c1, c2, c3) = self.getCell(mesh.vertex[t_in[i]])
                cells[i] = np.array([c1, c2, c3], dtype=np.int_)


                # create hash key h_i for v_i's cell
                hashes[i] = hash((c1, c2, c3))

            # if any v_i have same hash
            if ((hashes[0] == hashes[1]) or (hashes[0] == hashes[2]) or (hashes[1] == hashes[2])):
                # discard triangle entirely (continue)          
                continue

            
            for i in range(3):
                # for each vertex index v_i in triangle t_{in}
                v_i = t_in[i]
                
                if hashes[i] not in self.vertex_id:
                    # if hash h_i not found before:
                    if (self.next_vertex_id == 3890):
                        self.query_hash = hashes[i]
                    # add h_i to hash map vertex-id with key h_i and value "next-vertex-id"
                    self.vertex_id[hashes[i]] = self.next_vertex_id
                    self.next_vertex_id += 1

                    # create array in vertices-in-cell[h_i]
                    self.vertices_in_cell[hashes[i]] = []
                
                # add vertex index $v_i$ to array in vertices-in-cell[h_i]
                self.vertices_in_cell[hashes[i]].append(v_i)

            # initialize t_{out} as (vertex-id[h_1], vertex-id[h_2], vertex-id[h_3])
            t_out = np.zeros(3, dtype=np.int_)
            for i in range(3):
                t_out[i] = self.vertex_id[hashes[i]]
            
            sorted_t_out = np.sort(t_out)
            # add t_{out} to new mesh
            if (sorted_t_out[0], sorted_t_out[1], sorted_t_out[2]) in new_triangles_dict:
                # exact triangle already exists in new mesh
                continue
            
            new_triangles_dict[(sorted_t_out[0], sorted_t_out[1], sorted_t_out[2])] = 1
            new_triangles.append(t_out)
            face_colors.append(GREEN)
        
        simple_mesh = geomproc.mesh()
        
        simple_mesh.face = np.array(new_triangles)
        simple_mesh.fcolor = np.array(face_colors)

        # set new-vertices as the mesh's vertex array
        simple_mesh.vertex, approx_vnormals = self.calculateMergedVertices(mesh)
        simple_mesh.compute_vertex_and_face_normals()
        orientFaces(simple_mesh, approx_vnormals)

        return simple_mesh
        


def sq_norm(v):
    '''
    Compute squared norm of a 3D vector
    Source: geomproc library
    '''
    if (v.shape[0] == 2):
        return v[0]*v[0] + v[1]*v[1]
    else:
        return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]


def voronoiAngle(v1, v2, v3):
    '''
    Compute the area of the Voronoi region around v1 restricted to the face
    source: geomproc library
    '''
    
    vec1 = v1 - v2
    vec1 /= np.linalg.norm(vec1)
    vec2 = v1 - v3
    vec2 /= np.linalg.norm(vec2)
    cosine = np.dot(vec1, vec2)
    sine = np.linalg.norm(np.cross(vec1, vec2))
    if (sine == 0):
        print("WARN: triangle with sin(x)=0", v1, v2, v3)
    cot_v1 = cosine/sine
    if (cot_v1 < 0):
        # obtuse angle - return 1/2 area of triangle
        return np.linalg.norm(np.cross(v1 - v2, v1 - v3))/4

    
    # Compute the cotangents of the two corners opposite to v1
    vec1 = v3 - v2
    vec1 /= np.linalg.norm(vec1)
    vec2 = v1 - v2
    vec2 /= np.linalg.norm(vec2)
    cosine = np.dot(vec1, vec2)
    sine = np.linalg.norm(np.cross(vec1, vec2))
    cot_v2 = cosine/sine
    if (cot_v2 < 0):
        # obtuse angle - return 1/4 area of triangle
        return np.linalg.norm(np.cross(v3 - v2, v1 - v2))/8

    vec1 = v1 - v3
    vec1 /= np.linalg.norm(vec1)
    vec2 = v2 - v3
    vec2 /= np.linalg.norm(vec2)
    cosine = np.dot(vec1, vec2)
    sine = np.linalg.norm(np.cross(vec1, vec2))
    cot_v3 = cosine/sine
    if (cot_v3 < 0):
        # obtuse angle - return 1/4 area of triangle
        return np.linalg.norm(np.cross(v1 - v3, v2 - v3))/8

    # Compute the area based on cotangents and edge lengths

    # Voronoi area that touches (v1, v2)
    area_v2 = sq_norm(v1 - v2)*cot_v3
    # Voronoi area that touches (v1, v3)
    area_v3 = sq_norm(v1 - v3)*cot_v2
    # Total Voronoi area on this face for v1
    varea = (1/8)*(area_v2 + area_v3) 
    return varea

def flipTriangle(tri):
    '''
    Helper function for inverting a triangle's orientation
    '''
    tri[1], tri[0] = tri[0], tri[1]
    return tri

def calculateMeshVoronoi(mesh: geomproc.mesh):
    '''
    Calculate Voronoi area of all vertices in a mesh
    '''
    vareas = np.zeros(mesh.vertex.shape[0])

    for f in range(len(mesh.face)):
        for i in range(3):
            a = mesh.face[f][i]
            b = mesh.face[f][(i+1) % 3]
            c = mesh.face[f][(i+2) % 3]

            temp_varea = voronoiAngle(mesh.vertex[a], mesh.vertex[b], mesh.vertex[c])
            vareas[a] += temp_varea
    return vareas

def orientFaces(mesh : geomproc.mesh, approx_vnormals):
    '''
    Orient all triangle faces of a mesh so that they face the direction according to their vertices
    '''
    flipped = 0
    if (mesh.fnormal.shape[0] == 0):
        mesh.compute_vertex_and_face_normals()
    for f in range(mesh.fnormal.shape[0]):
        approx_fnormal = np.average(approx_vnormals[mesh.face[f]], axis=0)
        if (np.dot(mesh.fnormal[f], approx_fnormal) < 0):
            mesh.face[f] = flipTriangle(mesh.face[f])
            flipped += 1
    
    print(f"Flipped {flipped} triangles")

def gridSimplification(mesh: geomproc.mesh, num_cubes=1000, average_type=VORONOI_AVERAGE) -> geomproc.mesh:
    dimension = mesh.vertex.shape[1]
    if (len(mesh.vif) == 0):
        mesh.compute_connectivity()

    grid = Grid(num_cubes=num_cubes, dimension=dimension, average_type=average_type)
    simple_mesh = grid.simplifyMesh(mesh)

    return simple_mesh


def showGrid(num_cubes):
    '''
    Debugging method to save a mesh containing points, representing the corners of each cell in the Grid
    '''
    points = []
    grid = Grid(num_cubes=num_cubes)
    pc = geomproc.pcloud()
    step = 1/num_cubes
    for x in range(num_cubes*2 + 1):
        for y in range(num_cubes*2 + 1):
            for z in range(num_cubes*2 + 1):
                test_point = np.array([x,y,z]) * step
                test_point -= np.array([1,1,1])
                points.append(grid.getCellCorner(grid.getCell(test_point)))

    pc.point = np.array(points)
    wo = geomproc.write_options()
    wo.write_vertex_colors = True
    pts = geomproc.create_points(pc.point, radius=0.005, color=[0, 0, 1])
    # Combine everything together
    pc_mesh = geomproc.mesh()
    pc_mesh.append(pts)
    pc_mesh.save("./output/cube_corners.obj", wo)


if __name__ == "__main__":
    # AVERAGE_TYPES:
    # 0. VORONOI_AVERAGE
    # 1. MEAN_AVERAGE
    # 2. CORNER_AVERAGE

    average_type = VORONOI_AVERAGE
    num_cubes = 50
    # showGrid(num_cubes)

    # Specify input mesh
    # mesh_name = 'crater'
    # mesh_name = 'beetle'
    # mesh_name = 'bunny'
    # mesh_name = 'camel'
    # mesh_name = "cow"
    # mesh_name = "xyzrgb_dragon"
    # mesh_name = "kid"
    # mesh_name = "cat"
    mesh_name = "Armadillo"


    input_mesh = GEOMPROC + '/meshes/' + mesh_name + '.obj'
    output_mesh = 'output/' + mesh_name + f"_simplified_{AVERAGE_TYPES[average_type]}_{num_cubes}"
    output_mesh += ".ply"
    wo = geomproc.write_options()
    wo.write_face_colors = True
    mesh = geomproc.load(input_mesh)

    mesh.normalize()
    # Save normalized mesh
    mesh.save("./output/" + mesh_name + "_normalized.obj")
    mesh.name = mesh_name
    start_time = time.time()
    result = gridSimplification(mesh, num_cubes, average_type)
    result.save(output_mesh, wo)
    end_time = time.time()
    print('Execution time = ' + str(end_time - start_time) +'s')
