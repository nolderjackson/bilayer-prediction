import numpy as np
from numpy import sqrt, cos, sin, array

def angle(v1,v2):
    return np.rad2deg(np.arccos(np.dot(v1,v2)))
def VectorRotation(axis, angle):
    x, y, z = axis
    c = cos(angle)
    n = 1 - c
    s = sin(angle)
    rotation = array([[x*x*n + c, x*y*n - z*s, x*z*n + y*s],
                     [x*y*n + z*s, y*y*n + c, y*z*n - x*s],
                     [x*z*n - y*s, y*z*n + x*s, z*z*n + c]])
    return rotation

def MakeRandomVector():
    vec = np.random.rand(3)
    vec /= np.linalg.norm(vec)  # Normalize to unit length
    return vec

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def clean_vectors(vector_list):
    return [v for v in vector_list if v[2] >= 0]

def clean_pos_vectors(vector_list):
    vlist=[]
    for v in vector_list:
        if  v[2] >= 0:
            vlist.append(v)
        else:
            vlist.append(-1*v)
    return vlist

def calc_Vectors_In_Normal_Plane(vec, subdivisions):
    vector_list = []
    angle_list = np.deg2rad(np.arange((360/subdivisions), 180, (360/subdivisions))) #trust math works out generates 1/2 of all angles

    perp_vec = normalize(np.cross(vec, array([0,0,1])))
    if np.all(np.abs(perp_vec) < 1e-6): #if vec = norm then set perp = x-axis
        perp_vec = array([1,0,0])
    else:
        normalize(perp_vec)
    vector_list.extend((perp_vec, -perp_vec))
    for angle in angle_list:

        vectoapend = np.matmul(VectorRotation(vec, angle), perp_vec)
        vector_list.extend((vectoapend, -vectoapend))
    return vector_list

def MakeBindingVectors(dipole_vec, binding_angle, subdivisions):
    uvectors = calc_Vectors_In_Normal_Plane(dipole_vec, subdivisions)
    binding_vectors = [np.matmul(VectorRotation(vec, binding_angle), dipole_vec) for vec in uvectors]
    return binding_vectors

    #already includes linker radius
def ConstructBilayer(BottomSize:int, TopSize:int, BottomDistanceToDipole:int, TopDistanceToDipole:int,
                     BottomDipole, TopDipole, BottomBindingAngle:int, TopBindingAngle:int, Subdivisions:int):
    Bilayer_Index = []
    rlen_all = []
    rvec_all = []
    Top_Directions = []
    Bottom_Vec_all = []
    Top_Vec_all=[]
    i = 0
    Bottom_Directions = MakeBindingVectors(BottomDipole, BottomBindingAngle, Subdivisions)

    Bottom_Directions = clean_pos_vectors(Bottom_Directions)
    half_Top_Directions = MakeBindingVectors(TopDipole, TopBindingAngle, Subdivisions)
    negated_half = np.negative(half_Top_Directions)  # Returns a NumPy array
    Top_Directions.extend(half_Top_Directions)
    Top_Directions.extend(negated_half)

    BottomDistanceToDipole = [BottomDistanceToDipole*Bottom for Bottom in Bottom_Directions]
    TopDistanceToDipole = [TopDistanceToDipole*Top for Top in Top_Directions]

    Bottom_Sizes = [BottomSize*Bottom for Bottom in Bottom_Directions]
    Top_Sizes = [TopSize*Top for Top in Top_Directions]

    for Bindex in range(len(Bottom_Directions)):
        for Tindex in range(len(Top_Directions)):
            if Bottom_Sizes[Bindex][2] + Top_Sizes[Tindex][2] > Bottom_Sizes[Bindex][2]:
                Bilayer_Index.append([Bindex, Tindex])
            else:
                i += 1
    # print(Bilayer_Index)
     
    for pair in Bilayer_Index:
        rvec = BottomDistanceToDipole[pair[0]] + TopDistanceToDipole[pair[1]]
        Bottom_Vec_all.append(Bottom_Sizes[pair[0]])
        Top_Vec_all.append(Top_Sizes[pair[1]])
        rlen_all.append(np.linalg.norm(rvec))
        rvec_all.append(normalize(rvec))


    return rvec_all, rlen_all, Bottom_Vec_all, Top_Vec_all

def calcmuD(theta_nmuD):
    return array([sin(theta_nmuD), 0, cos(theta_nmuD)])

def calcmuA(theta_beta, theta_nmuD, theta_nmuA):
    if theta_nmuA == 0:
        return array([0,0,1])
    elif theta_nmuD == 0:
        return array([sin(theta_nmuD), 0, cos(theta_nmuD)])
    else: 
        try:
            return array([
                (cos(theta_beta) - cos(theta_nmuA) * cos(theta_nmuD)) / sin(theta_nmuD),
                -sqrt(-((cos(theta_beta) - cos(theta_nmuA) * cos(theta_nmuD)) / (sin(theta_nmuA) * sin(theta_nmuD)))**2 + 1) * sin(theta_nmuA),
                cos(theta_nmuA)])
        except:
            theta_nmuA = np.round(theta_nmuA, 3)
            theta_beta = np.round(theta_beta, 3)
            if theta_nmuA == 0:
                return array([0,0,1])
            elif theta_nmuD == 0:
                return array([sin(theta_nmuD), 0, cos(theta_nmuD)])
        
def calcKappa2(rvec_list, theta_beta, muA, muD):
    kappa2_list = []
    Cos_Beta = cos(theta_beta)
    for rvec in rvec_list:
        muDdot= np.dot(muD, rvec)
        muAdot= np.dot(muA, rvec)
        if muDdot < 0:
            muDdot = -muDdot
        if muDdot < 0:
            muAdot = -muAdot
        kappa2 = (Cos_Beta - (3*np.dot(muD, rvec)*np.dot(muA, rvec)))**2
        kappa2_list.append(kappa2)
    return kappa2_list



import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def visualize_bilayer_vectors(dipoleD, dipoleA, BottomBindingAngle, TopBindingAngle, Subdivisions,
                              BottomSize, TopSize, BottomDistanceToDipole, TopDistanceToDipole,
                              dipole_scale=1):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Vector Stack Visualization with Dipole Directions")
    Top_Vecs = []
    # Generate and clean binding vectors
    Bottom_Vecs = MakeBindingVectors(dipoleD, BottomBindingAngle, Subdivisions)
    half_Top_Vecs = MakeBindingVectors(dipoleA, TopBindingAngle, Subdivisions)
    Bottom_Vecs = clean_vectors(Bottom_Vecs)

    negated_half = np.negative(half_Top_Vecs)  # Returns a NumPy array
    Top_Vecs.extend(half_Top_Vecs)
    Top_Vecs.extend(negated_half)
    pair_count = 0

    for i, bvec in enumerate(Bottom_Vecs):
        for j, tvec in enumerate(Top_Vecs):
            # if BottomSize * bvec[2] + TopSize * tvec[2] <= 0:
            #     continue

            bvec = normalize(bvec)
            tvec = normalize(tvec)

            origin = np.array([0, 0, 0])
            joint_point = origin + BottomSize * bvec
            dipoleD_pos = origin + (BottomSize - BottomDistanceToDipole) * bvec
            dipoleA_pos = joint_point + TopDistanceToDipole * tvec

            # Bottom binding vector
            ax.quiver(*origin, *(BottomSize * bvec), color='blue', label='Bottom Vector' if pair_count == 0 else "")

            # Top binding vector
            ax.quiver(*joint_point, *(TopSize * tvec), color='red', label='Top Vector' if pair_count == 0 else "")

            # Dipole-D and Dipole-A markers
            ax.scatter(*dipoleD_pos, color='blue', marker='o', label='Dipole D' if pair_count == 0 else "")
            ax.scatter(*dipoleA_pos, color='red', marker='^', label='Dipole A' if pair_count == 0 else "")

            # Dipole direction vectors
            ax.quiver(*dipoleD_pos, *(dipoleD * dipole_scale), color='cyan', label='DipoleD dir' if pair_count == 0 else "")
            ax.quiver(*dipoleA_pos, *(dipoleA * dipole_scale), color='magenta', label='DipoleA dir' if pair_count == 0 else "")

            # Dipole-Dipole link
            ax.plot([dipoleD_pos[0], dipoleA_pos[0]],
                    [dipoleD_pos[1], dipoleA_pos[1]],
                    [dipoleD_pos[2], dipoleA_pos[2]],
                    color='green', linestyle='--', linewidth=1, label='Dipole-Dipole Link' if pair_count == 0 else "")

            pair_count += 1

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_box_aspect([1, 1, 1])
    ax.legend()
    plt.show()


if __name__ == '__main__':

    # Define angles----------------------------------------------------------------------------------
    theta_beta = np.deg2rad(50)
    theta_nmuBL = np.deg2rad(31)
    theta_nmuTL = np.deg2rad(50)
    theta_AmuA = np.deg2rad(90)
    theta_DmuD = np.deg2rad(75)
    BottomDistanceToDipole=1
    TopDistanceToDipole=1.5
    LengthBL = 1   # Arb Unit
    LengthTL = 1
    radius = 0
    subdivisions = 10 #Smallest num to get 1,000,000 total structures

    """
    subdivisions O(n^3)
    """
    #---------------------------------------------------------------------------------------------------

    LengthTL = LengthTL + radius
    LengthBL = LengthBL + radius


    muBL = calcmuD(theta_nmuBL)
    muTL = calcmuA(theta_beta, theta_nmuBL, theta_nmuTL)

    vec_all, rlen_all, bvecall, tvecall = ConstructBilayer(LengthBL, LengthTL, BottomDistanceToDipole, TopDistanceToDipole, muBL, muTL, theta_DmuD, theta_AmuA, subdivisions)
    # visualize_bilayer_vectors(dipoleD, dipoleA, BottomBindingAngle, TopBindingAngle, Subdivisions,
    #                           BottomSize, TopSize, BottomDistanceToDipole, TopDistanceToDipole,
    #                           dipole_scale=1)
    for t in tvecall:
        print(angle(t, muTL))
    for b in bvecall:
        print(angle(b, muBL))    
