import numpy as np
from numpy import cos, sin, deg2rad

def MakeAngles(theta_nmuBL,
               Beta_min, Beta_Max, Beta_Subdivisions,
               theta_nmuTL_min, theta_nmuTL_max, theta_nmuTL_Subdivisions):
    print(theta_nmuBL,
               Beta_min, Beta_Max, Beta_Subdivisions,
               theta_nmuTL_min, theta_nmuTL_max, theta_nmuTL_Subdivisions)
    if Beta_min < 0 or theta_nmuTL_min < 0:
        return
    if Beta_Max > 90 or theta_nmuTL_max > 90:
        return
    if theta_nmuBL < 0 or theta_nmuBL > 90:
        return
    
    AngleList = []
    BadAngleList = []

    Beta_range = np.linspace(Beta_min, Beta_Max, Beta_Subdivisions)
    theta_nmuA_range = np.linspace(theta_nmuTL_min, theta_nmuTL_max, theta_nmuTL_Subdivisions)

    theta_nmuDr = np.deg2rad(theta_nmuBL)
    cosnmuD = np.cos(theta_nmuDr)
    sin2nmuD = np.sin(theta_nmuDr) ** 2

    for beta in Beta_range:

        betar = deg2rad(beta)

        cosB = cos(betar)

        for theta_nmuA in theta_nmuA_range:

            theta_nmuAr = deg2rad(theta_nmuA)

            cosnmuA = cos(theta_nmuAr)
            sin2nmuA = sin(theta_nmuAr)**2

            if theta_nmuBL == 0 and beta == theta_nmuA:
                AngleList.append(f"{beta:.1f},{theta_nmuBL:.1f},{theta_nmuA:.1f},{0}")
            elif theta_nmuBL == 0 and beta != theta_nmuA:
                BadAngleList.append(f"{beta:.1f},{theta_nmuBL:.1f},{theta_nmuA:.1f},{0}")

            if theta_nmuA == 0 and beta == theta_nmuBL:
                AngleList.append(f"{beta:.1f},{theta_nmuBL:.1f},{theta_nmuA:.1f},{0}")
            elif theta_nmuA == 0 and beta != theta_nmuBL:
                BadAngleList.append(f"{beta:.1f},{theta_nmuBL:.1f},{theta_nmuA:.1f},{0}")

            if sin2nmuD * sin2nmuA != 0:
                numerator = (cosB - (cosnmuD * cosnmuA)) ** 2
                denominator = sin2nmuD * sin2nmuA
                fraction = numerator / denominator
                if fraction <= 1:
                    val = np.sqrt(1 - fraction)
                    AngleList.append(f"{beta:.1f},{theta_nmuBL:.1f},{theta_nmuA:.1f},{val:.5f}")
                else:
                    BadAngleList.append(f"{beta:.1f},{theta_nmuBL:.1f},{theta_nmuA:.1f},{fraction:.5f}")

    with open("AngleDict.csv", "w") as f:
        f.write("beta,theta_nmuBL,theta_nmuTL,value\n")
        f.write("\n".join(AngleList))

    with open("BadAngleDict.csv", "w") as f:
        f.write("beta,theta_nmuBL,theta_nmuTL,value\n")
        f.write("\n".join(BadAngleList))
    return

if __name__ == "__main__":

    theta_nmuBL = 31
    
    beta_min = 0
    beta_max = 90
    beta_subdivisions = 19

    theta_nmuTL_min = 0
    theta_nmuTL_max = 90
    theta_nmuTL_subdivisions = 19
    MakeAngles(theta_nmuBL,
               beta_min, beta_max, beta_subdivisions,
               theta_nmuTL_min, theta_nmuTL_max, theta_nmuTL_subdivisions)

    # O(n^2)
    # if n = num of subdivisions
    # recomended range 9-19
    # min to max is 0 indexed