import cmb_ilc
reload(cmb_ilc)
from cmb_ilc import *


##################################################################################

# Specifications
Nu = np.array([90.e9, 150.e9])   # [Hz]
Beam = np.array([1.4, 1.4])   # [arcmin]
Noise = np.array([14., 20.])  # [muK*arcmin]



##################################################################################
# init

cmbIlc = CMBILC(Nu, Beam, Noise)

cmbIlc.plotWeightsIlcCmb()
cmbIlc.plotPowerIlcCmb()

cmbIlc.plotWeightsIlcTsz()
cmbIlc.plotPowerIlcTsz()

cmbIlc.plotWeightsConstrainedIlcTszNoCib()
cmbIlc.plotPowerConstrainedIlcTszNoCib()
