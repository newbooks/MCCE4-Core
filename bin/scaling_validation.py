#!/usr/bin/env python
"""
MCCE4 Tool: Validate Scaling Factors
Coulomb potential is propotional to the charge, we need to valid PB potential would be propotional to the charge as well.
If this is the case, we only need to scale the Coulomb potential by a factor to get the PB potential.
"""

PDB_Template = """
ATOM   1909  N   ARG A0128_000  -3.675   8.073 -18.684   1.550       0.000      BKO000_000
ATOM   1910  CA  ARG A0128_000  -5.057   7.918 -19.082   1.700       0.000      BKO000_000
ATOM   1911  C   ARG A0128_000  -5.840   7.425 -17.887   1.700       0.000      BKO000_000
ATOM   1912  O   ARG A0128_000  -5.606   6.382 -17.247   1.520       0.000      BKO000_000
ATOM   1913  H   ARG A0128_000  -2.980   7.281 -18.963   1.200       0.000      BKO000_000
ATOM   1914  HA  ARG A0128_000  -5.413   8.881 -19.448   1.200       0.000      BKO000_000
ATOM   1915  CB  ARG A0128_001  -5.216   6.873 -20.190   1.700       0.000      01G000_000
ATOM   1916  CG  ARG A0128_001  -4.052   6.385 -21.007   1.700       0.000      01G000_000
ATOM   1917  CD  ARG A0128_001  -4.343   5.157 -21.845   1.700       0.000      01G000_000
ATOM   1918  NE  ARG A0128_001  -4.390   5.271 -23.247   1.550       0.000      01G000_000
ATOM   1919  CZ  ARG A0128_001  -4.334   5.521 -24.503   1.700       0.000      01G000_000
ATOM   1920  NH1 ARG A0128_001  -3.778   6.614 -25.049   1.550       0.000      01G000_000
ATOM   1921  NH2 ARG A0128_001  -4.874   4.635 -25.349   1.550       0.000      01G000_000
ATOM   1922  HB2 ARG A0128_001  -5.932   7.289 -20.899   1.200       0.000      01G000_000
ATOM   1923  HB3 ARG A0128_001  -5.645   5.990 -19.716   1.200       0.000      01G000_000
ATOM   1924  HG2 ARG A0128_001  -3.235   6.147 -20.326   1.200       0.000      01G000_000
ATOM   1925  HG3 ARG A0128_001  -3.745   7.189 -21.676   1.200       0.000      01G000_000
ATOM   1926  HD2 ARG A0128_001  -5.313   4.777 -21.525   1.200       0.000      01G000_000
ATOM   1927  HD3 ARG A0128_001  -3.571   4.424 -21.612   1.200       0.000      01G000_000
ATOM   1928 HH11 ARG A0128_001  -3.782   6.745 -26.131   1.200       0.000      01G000_000
ATOM   1929 HH12 ARG A0128_001  -3.323   7.371 -24.410   1.200       0.000      01G000_000
ATOM   1930 HH21 ARG A0128_001  -4.848   4.821 -26.423   1.200       0.000      01G000_000
ATOM   1931 HH22 ARG A0128_001  -5.340   3.728 -24.964   1.200       0.000      01G000_000
"""

CHARGE_SOURCE = " CB "      # The source of charge for the potential calculation, will change in range [-1, 1]
POTENTIAL_SITE = " NE "     # The site of potential to be reported, maintain constant charge 1

