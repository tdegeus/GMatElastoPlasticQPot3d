
import GMatElastoPlasticQPot3d.Cartesian3d as GMat
import numpy as np

# ==================================================================================================

def EQ(a,b):
  assert np.abs(a-b) < 1.e-12

# ==================================================================================================

# material model
# - parameters
kappa = 12.3
mu    = 45.6
# - model
mat = GMat.Elastic(kappa,mu)

# simple shear + volumetric deformation
# - parameters
gamma = 0.02
epsm  = 0.12
# - strain
Eps = [[ epsm  , gamma, 0.  ],
       [ gamma , epsm , 0.  ],
       [ 0.    , 0.   , epsm]]
# - stress
Sig = mat.Sig(Eps)
# - analytical solution
EQ( Sig[0,0], 3. * kappa * epsm  )
EQ( Sig[1,1], 3. * kappa * epsm  )
EQ( Sig[2,2], 3. * kappa * epsm  )
EQ( Sig[0,1], 2. * mu    * gamma )
EQ( Sig[1,0], 2. * mu    * gamma )
EQ( Sig[0,2], 0.                 )
EQ( Sig[2,0], 0.                 )
EQ( Sig[1,2], 0.                 )
EQ( Sig[2,1], 0.                 )
# - plastic strain
EQ( mat.epsp(Eps), 0 )
# - yield strain index
EQ( mat.find(Eps), 0 )

# ==================================================================================================

# material model
# - parameters
kappa = 12.3
mu    = 45.6
# - model
mat = GMat.Cusp(kappa,mu,[0.01, 0.03, 0.10])

# simple shear + volumetric deformation
# - parameters
gamma = 0.02
epsm  = 0.12
# - strain
Eps = [[ epsm  , gamma, 0.  ],
       [ gamma , epsm , 0.  ],
       [ 0.    , 0.   , epsm]]
# - stress
Sig = mat.Sig(Eps)
# - analytical solution
EQ( Sig[0,0], 3. * kappa * epsm  )
EQ( Sig[1,1], 3. * kappa * epsm  )
EQ( Sig[2,2], 3. * kappa * epsm  )
EQ( Sig[0,1], 0.                 )
EQ( Sig[1,0], 0.                 )
EQ( Sig[0,2], 0.                 )
EQ( Sig[2,0], 0.                 )
EQ( Sig[1,2], 0.                 )
EQ( Sig[2,1], 0.                 )
# - plastic strain
EQ( mat.epsp(Eps), 0.02 )
# - yield strain index
EQ( mat.find(Eps), 1 )

# ==================================================================================================

# material model
# - parameters
kappa = 12.3
mu    = 45.6
# - model
mat = GMat.Smooth(kappa,mu,[0.01, 0.03, 0.10])

# simple shear + volumetric deformation
# - parameters
gamma = 0.02
epsm  = 0.12
# - strain
Eps = [[ epsm  , gamma, 0.  ],
       [ gamma , epsm , 0.  ],
       [ 0.    , 0.   , epsm]]
# - stress
Sig = mat.Sig(Eps)
# - analytical solution
EQ( Sig[0,0], 3. * kappa * epsm  )
EQ( Sig[1,1], 3. * kappa * epsm  )
EQ( Sig[2,2], 3. * kappa * epsm  )
EQ( Sig[0,1], 0.                 )
EQ( Sig[1,0], 0.                 )
EQ( Sig[0,2], 0.                 )
EQ( Sig[2,0], 0.                 )
EQ( Sig[1,2], 0.                 )
EQ( Sig[2,1], 0.                 )
# - plastic strain
EQ( mat.epsp(Eps), 0.02 )
# - yield strain index
EQ( mat.find(Eps), 1 )

# ==================================================================================================

# parameters
kappa = 12.3
mu    = 45.6

# allocate matrix
mat = GMat.Matrix([3,2])

# row 0: elastic
I      = np.zeros(mat.shape(), dtype='int')
I[0,:] = 1
mat.setElastic(I,kappa,mu)

# row 1: cups
I      = np.zeros(mat.shape(), dtype='int')
I[1,:] = 1
mat.setCusp(I,kappa,mu,[0.01,0.03,0.10])

# row 2: smooth
I      = np.zeros(mat.shape(), dtype='int')
I[2,:] = 1
mat.setSmooth(I,kappa,mu,[0.01,0.03,0.10])

# simple shear + volumetric deformation
# - parameters
gamma = 0.02;
epsm  = 0.12;
# - strain
Eps          = np.zeros((3,2,3,3))
Eps[:,:,0,0] = epsm
Eps[:,:,1,1] = epsm
Eps[:,:,2,2] = epsm
Eps[:,:,0,1] = gamma
Eps[:,:,1,0] = gamma
# - stress & plastic strain
Sig  = mat.Sig (Eps)
epsp = mat.epsp(Eps)

# - analytical solution
EQ( Sig[0,0,0,0], 3. * kappa * epsm ); EQ( Sig[0,1,0,0], 3. * kappa * epsm )
EQ( Sig[0,0,1,1], 3. * kappa * epsm ); EQ( Sig[0,1,1,1], 3. * kappa * epsm )
EQ( Sig[0,0,0,1], 2. * mu    * gamma); EQ( Sig[0,1,0,1], 2. * mu    * gamma)
EQ( Sig[0,0,0,1], 2. * mu    * gamma); EQ( Sig[0,1,1,0], 2. * mu    * gamma)
EQ( Sig[0,0,0,2], 0.                ); EQ( Sig[0,1,0,2], 0.                )
EQ( Sig[0,0,2,0], 0.                ); EQ( Sig[0,1,2,0], 0.                )
EQ( Sig[0,0,1,2], 0.                ); EQ( Sig[0,1,1,2], 0.                )
EQ( Sig[0,0,2,1], 0.                ); EQ( Sig[0,1,2,1], 0.                )
# -
EQ( Sig[1,0,0,0], 3. * kappa * epsm ); EQ( Sig[1,1,0,0], 3. * kappa * epsm )
EQ( Sig[1,0,1,1], 3. * kappa * epsm ); EQ( Sig[1,1,1,1], 3. * kappa * epsm )
EQ( Sig[1,0,0,1], 0.                ); EQ( Sig[1,1,0,1], 0.                )
EQ( Sig[1,0,0,1], 0.                ); EQ( Sig[1,1,1,0], 0.                )
EQ( Sig[1,0,0,2], 0.                ); EQ( Sig[1,1,0,2], 0.                )
EQ( Sig[1,0,2,0], 0.                ); EQ( Sig[1,1,2,0], 0.                )
EQ( Sig[1,0,1,2], 0.                ); EQ( Sig[1,1,1,2], 0.                )
EQ( Sig[1,0,2,1], 0.                ); EQ( Sig[1,1,2,1], 0.                )
# -
EQ( Sig[2,0,0,0], 3. * kappa * epsm ); EQ( Sig[2,1,0,0], 3. * kappa * epsm )
EQ( Sig[2,0,1,1], 3. * kappa * epsm ); EQ( Sig[2,1,1,1], 3. * kappa * epsm )
EQ( Sig[2,0,0,1], 0.                ); EQ( Sig[2,1,0,1], 0.                )
EQ( Sig[2,0,0,1], 0.                ); EQ( Sig[2,1,1,0], 0.                )
EQ( Sig[2,0,0,2], 0.                ); EQ( Sig[2,1,0,2], 0.                )
EQ( Sig[2,0,2,0], 0.                ); EQ( Sig[2,1,2,0], 0.                )
EQ( Sig[2,0,1,2], 0.                ); EQ( Sig[2,1,1,2], 0.                )
EQ( Sig[2,0,2,1], 0.                ); EQ( Sig[2,1,2,1], 0.                )
# - plastic strain
EQ( epsp[0,0], 0    ); EQ( epsp[0,1], 0    )
EQ( epsp[1,0], gamma); EQ( epsp[1,1], gamma)
EQ( epsp[2,0], gamma); EQ( epsp[2,1], gamma)

# ==================================================================================================

print('All checks passed')
