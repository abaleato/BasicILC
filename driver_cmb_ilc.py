import cmb_ilc
reload(cmb_ilc)
from cmb_ilc import *


##################################################################################

# Specifications
Nu = np.array([27.e9, 39.e9, 93.e9, 145.e9, 225.e9, 280.e9])   # [Hz]
Beam = np.array([7.4, 5.1, 2.2, 1.4, 1.0, 0.9])   # [arcmin]
Noise = np.array([52., 27., 5.8, 6.3, 15., 37.])  # [muK*arcmin]



##################################################################################
# init

cmbIlc = CMBILC(Nu, Beam, Noise)

cmbIlc.plotWeightsIlcCmb()
cmbIlc.plotPowerIlcCmb()

cmbIlc.plotWeightsIlcTsz()
cmbIlc.plotPowerIlcTsz()

cmbIlc.plotWeightsConstrainedIlcTszNoCib()
cmbIlc.plotPowerConstrainedIlcTszNoCib()
