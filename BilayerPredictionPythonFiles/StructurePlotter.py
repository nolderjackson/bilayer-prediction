import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
def StructurePlot(beta, theta_nmuTL, linkersize, filename):
    #Water2BilayerData.csv
    drop_cols = [
        'value',
        'Reorganization_Energy_avg',
        'G_eT_avg',
        'Attenuated_Coupling_avg',
        'kappa2_avg',
        'r_length_avg'
    ]

    base_path = Path(__file__).resolve().parent

    while True:
        data_path = base_path / (filename)
        try:
            df = pd.read_csv(data_path).drop(columns=drop_cols)
            break  # Success
        except FileNotFoundError:
            print("File not found. Please try again.\n")
        except pd.errors.EmptyDataError:
            print("File is empty. Please try again.\n")
        except PermissionError:
            print("No file named. Please try again.\n")

    # get a list of the columns
    cols = list(df.columns)
    i, j = cols.index("theta_nmuD"), cols.index("theta_nmuA")
    cols[i], cols[j] = cols[j], cols[i]  # swap positions
    df = df[cols]

    print(df, '\n')

    def undodot(v1,v2):
        return np.rad2deg(np.arccos(np.dot(v1,v2)))
    
    def set_axes_equal_with_padding(ax, pad=1.0):
        """Set 3D plot axes to equal scale with optional padding."""
        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # use the largest range so scaling is equal
        plot_radius = 0.5 * max([x_range, y_range, z_range]) + pad

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    def calcmuBL(theta_nmuD):
        theta_nmuD = np.deg2rad(theta_nmuD)
        return np.array([np.sin(theta_nmuD), 0, np.cos(theta_nmuD)])

    def calcmuTL(theta_beta, theta_nmuD, theta_nmuA):
        theta_beta, theta_nmuD, theta_nmuA = np.deg2rad(theta_beta), np.deg2rad(theta_nmuD), np.deg2rad(theta_nmuA)
        if theta_nmuA == 0:
            return np.array([0,0,1])
        elif theta_nmuD == 0:
            return np.array([np.sin(theta_nmuD), 0, np.cos(theta_nmuD)])
        else: 
            try:
                return np.array([
                    (np.cos(theta_beta) - np.cos(theta_nmuA) * np.cos(theta_nmuD)) / np.sin(theta_nmuD),
                    -np.sqrt(-((np.cos(theta_beta) - np.cos(theta_nmuA) * np.cos(theta_nmuD)) / (np.sin(theta_nmuA) * np.sin(theta_nmuD)))**2 + 1) * np.sin(theta_nmuA),
                    np.cos(theta_nmuA)])
            except:
                theta_nmuA = np.round(theta_nmuA, 3)
                theta_beta = np.round(theta_beta, 3)
                if theta_nmuA == 0:
                    return np.array([0,0,1])
                elif theta_nmuD == 0:
                    return np.array([np.sin(theta_nmuD), 0, np.cos(theta_nmuD)])
                
    def normalize(v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm

    def CalcDipLoc(v, lsize):
        vnorm = normalize(v)
        vtrue = v - (vnorm*lsize)
        return vtrue/2

    def AssignVariables(beta_in, theta_nmuA_in):

        idx = ((df['beta'] - beta_in).abs() + (df['theta_nmuA'] - theta_nmuA_in).abs()).idxmin()

        # Get the row
        row = df.loc[idx]

        if row.empty:
            print('\nRow found to be empty, this error should never appear\n')
            row = None
 
        print(row,'\n')
        return row

    def CalcVecs(row):
        
        Bottom_FRET = np.fromstring(((row['Bottom_FRET'])).strip("[]"), sep=" ")
        Top_FRET = np.fromstring(((row['Top_FRET'])).strip("[]"), sep=" ")
        Bottom_eT = np.fromstring(((row['Bottom_eT'])).strip("[]"), sep=" ")
        Top_eT = np.fromstring(((row['Top_eT'])).strip("[]"), sep=" ")

        Vec1 = np.array(Bottom_FRET)
        Vec2 = np.array(Top_FRET)
        Vec3 = np.array(Bottom_eT)
        Vec4 = np.array(Top_eT)

        return Vec1, Vec2, Vec3, Vec4

    def CalcDips(row):
        beta = float((row['beta']))
        theta_nmuBL = float((row['theta_nmuD']))
        theta_nmuTL = float((row['theta_nmuA']))

        muBL = calcmuBL(theta_nmuBL)
        muTL = calcmuTL(beta, theta_nmuBL, theta_nmuTL)
        return muBL, muTL

    def PlotVectors(v1, v2, v3, v4, d1, d2, lsize, row):
        fig = plt.figure()
        origin = [0,0,0]
        plt.rcParams['savefig.dpi'] = 300


        plt.rcParams['font.size'] = 12  # Default font size for all text
        # plt.rcParams['axes.titlesize'] = 14  # Font size for plot titles
        plt.rcParams['axes.labelsize'] = 22  # Font size for x and y axis labels
        plt.rcParams['xtick.labelsize'] = 14  # Font size for x-axis tick labels
        plt.rcParams['ytick.labelsize'] = 14  # Font size for y-axis tick labels
        plt.rcParams['legend.fontsize'] = 20  # Font size for legend text
        plt.rcParams['figure.titlesize'] = 64 # Font size for figure titles
        d1locF = CalcDipLoc(v1, lsize)
        d2locF = CalcDipLoc(v2, lsize)
        d1locE = CalcDipLoc(v3, lsize)
        d2locE = CalcDipLoc(v4, lsize)

        fig.text(0.02,0.9,r'$\beta =$' + rf'{row['beta']:.1f}' + r'$\degree$'
                            '\n' + r'$\theta_{n, \mu_{TL}} =$' + rf'{row['theta_nmuA']:.1f}' + r'$\degree$', fontsize = 14)

        ax1 = fig.add_subplot(121, projection='3d')
        ax1.set_title(r"$\mathrm{Max}: \log{k_{FRET}} =$" + rf'{np.log10((row['kFRET_max'])):.2f}', pad=15, fontsize=14)
        ax1.text2D(0.5, -0.1, r"$\angle \hat{r}_{BL},\hat{r}_{TL} =$" + f'{undodot(normalize(-v1), normalize(v2)):.1f}°', 
            transform=ax1.transAxes,
            ha='center', va='top', fontsize = 14)
        ax1.quiver(*origin, *v1)
        ax1.quiver(*v1, *v2)

        ax1.quiver(*d1locF, *(2*d1), color='red')
        ax1.quiver(*d1locF, *(-2*d1), color='red')
        ax1.quiver(*(v1+v2-d2locF), *(2*d2), color='red')
        ax1.quiver(*(v1+v2-d2locF), *(-2*d2), color='red')

        ax2 = fig.add_subplot(122, projection='3d')
        ax2.set_title(r"$\mathrm{Max}: \log{k_{eT}} =$" + f'{np.log10((row['keT_max'])):.2f}', pad=15, fontsize=14)
        ax2.text2D(0.5, -0.1, r"$\angle \hat{r}_{BL},\hat{r}_{TL} =$" + f'{undodot(normalize(-v3), normalize(v4)):.1f}°', 
            transform=ax2.transAxes,
            ha='center', va='top', fontsize = 14)
        ax2.quiver(*origin, *v3)
        ax2.quiver(*v3, *v4)

        ax2.quiver(*d1locE, *(2*d1), color='red')
        ax2.quiver(*d1locE, *(-2*d1), color='red')
        ax2.quiver(*(v3+v4-d2locE), *(2*d2), color='red')
        ax2.quiver(*(v3+v4-d2locE), *(-2*d2), color='red')
        
        set_axes_equal_with_padding(ax1, pad=1.0)
        set_axes_equal_with_padding(ax2, pad=1.0)

        print('Angle between molecules for maximum', 'k_(FRET)', f'{min(undodot(normalize(v1), normalize(v2)), (undodot(normalize(-v2), normalize(v1))))}', '°')
        print('Angle between molecules for maximum', 'k_(eT)', f'{min(undodot(normalize(v3), normalize(v4)), (undodot(normalize(-v4), normalize(v3))))}', '°', '\n')

        plt.show()

    def RunProgram(beta, theta_nmuA):
        
        row = AssignVariables(beta_in=beta, theta_nmuA_in=theta_nmuA)
        vecs = CalcVecs(row)
        dips = CalcDips(row)
        # print(f'{dips=}')
        print("Plotting data\n")
        PlotVectors(*vecs, *dips, linkersize, row)

    RunProgram(beta, theta_nmuTL)


if __name__ == "__main__":

    filename = "BilayerData.csv"
    beta = 60
    theta_nmuTL = 30
    linker_size = 2              # in angstrom
    StructurePlot(beta, theta_nmuTL, linker_size, filename) 
    