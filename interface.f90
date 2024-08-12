MODULE interface

    IMPLICIT NONE
		!------------------------------------- LATTICE PARAMETERS:
		INTEGER(8), PARAMETER :: DIMHILB=2 !dimensão do espaço de Hilbert
		INTEGER(8), DIMENSION(DIMHILB) :: b1=(/ 1,0 /) !vetor da rede
		INTEGER(8), DIMENSION(DIMHILB) :: b2=(/ 0,1 /) !vetor da rede
		REAL(8), DIMENSION(DIMHILB) :: veck !vetor k
		INTEGER(8) :: itemp
		
		!------------------------------------- CONSTANTS AND INDEX PARAMETERS:
		REAL(8), PARAMETER :: zero=0.d0
		REAL(8), PARAMETER :: pi=DACOS(-1.d0)
		REAL(8), PARAMETER :: twopi=2.d0*DACOS(-1.d0)
		REAL(8), PARAMETER :: t=1.d0 !tight-binding parameter
		!REAL(8), PARAMETER :: U=1/0.077 !temporario para teste
		!REAL(8), PARAMETER :: ne=1.6 !temporario para teste
		REAL(8) :: delta !delta use in seed
		INTEGER(8), PARAMETER :: Nx=50  ! numero de vetores k em x - >256
		INTEGER(8), PARAMETER :: Ny=50  ! numero de vetores k em y - >256
		
		INTEGER(8), PARAMETER :: edgeinix=-Nx !ZB x-axis
		INTEGER(8), PARAMETER :: edgefimx=Nx !ZB x-axis
		
		INTEGER(8), PARAMETER :: edgeiniy=-Ny !ZB y-axis
		INTEGER(8), PARAMETER :: edgefimy=Ny !ZB y-axis
		
		INTEGER(8), PARAMETER :: Nk = (2*Nx+1)*(2*Ny+1)-(2*Nx+1)-(2*Ny+1)+1 !Para normalização eliminando bordas de cima e da direita
		INTEGER(8), PARAMETER :: edgebordax=Nx-1
		INTEGER(8), PARAMETER :: edgeborday=Ny-1
		COMPLEX(8), PARAMETER :: i=DCMPLX(0.d0,1.d0) ! numero complexo
		INTEGER(8) :: idimhilb ! indice que corre o espaço de Hilbert
		REAL(8):: BETA       !       INVERSE TEMPERATURE k_B*T 
		INTEGER :: initial=1
		REAL(8), PARAMETER :: kb=0.00008617d0 !Boltzmann constant in eV/K
		REAL(8), PARAMETER :: temperature=300.d0!10**(-4.D0)   !Kelvin
		REAL(8), DIMENSION(0:(4*Ny),DIMHILB) :: vecpath !vetor k
		
		!------------------------------------- HAMILTONIAN PARAMETERS:
		COMPLEX(8), DIMENSION(DIMHILB,DIMHILB) :: hdw
		COMPLEX(8), DIMENSION(DIMHILB,DIMHILB) :: hup 
		REAL(8), DIMENSION(edgeinix:edgefimx,edgeiniy:edgefimy,DIMHILB) :: edw
		REAL(8), DIMENSION(edgeinix:edgefimx,edgeiniy:edgefimy,DIMHILB) :: eup
		COMPLEX(8), DIMENSION(edgeinix:edgefimx,edgeiniy:edgefimy,DIMHILB,DIMHILB) :: vdw
		COMPLEX(8), DIMENSION(edgeinix:edgefimx,edgeiniy:edgefimy,DIMHILB,DIMHILB) :: vup
		COMPLEX(8) :: gama
		
		REAL(8), DIMENSION(0:(4*Nx),DIMHILB) :: edw_path
		REAL(8), DIMENSION(0:(4*Nx),DIMHILB) :: eup_path
	
		!------------------------------------- ENERGY PARAMETERS:
		REAL(8), DIMENSION(DIMHILB) :: minenergy ! minimum value of eigenenergy edw and eup
		REAL(8), DIMENSION(DIMHILB) :: maxenergy ! maximum value of eigenenergy edw and eup
		INTEGER(8), PARAMETER :: Nenergy=500*(Nk)! numero de pontos na energy
		INTEGER(8) :: ienergy ! indice que corre a energia
		REAL(8) :: step_energy ! passo da energy
		real(8), DIMENSION(0:Nenergy) :: w  !range de energy
		INTEGER(8) :: nsteps
		REAL(8), PARAMETER :: val_step_energy=0.001d0
		
		!------------------------------------- U AND NE PARAMETERS:
		INTEGER(8), PARAMETER :: NVALtU=15 !Quantos valores estamos calculando de U
		INTEGER(8) :: itu
		REAL(8) :: tU !Variable t/U - Hubbard
		REAL(8) :: U !Variable U - Hubbard
		REAL(8),PARAMETER :: tUMAX=0.2D0 !Valor max do U
		INTEGER(8), DIMENSION(0:NVALtU) :: NVECtU 
		REAL(8),PARAMETER :: tU0 = 0.1d0 !Valor min de U
		REAL(8) :: steptU
		
		INTEGER(8), PARAMETER :: NVALne=15 !Quantos valores estamos calculando de ne
		INTEGER(8) :: ine
		REAL(8) :: ne !Number of electrons 
		REAL(8),PARAMETER :: NEMAX=1.999d0 !Valor max de ne
		REAL(8),PARAMETER :: ne0 = 0.2d0 !Valor min de ne
		INTEGER(8), DIMENSION(0:NVALne) :: NVECne 
		REAL(8) :: stepne
		REAL(8), DIMENSION(0:Nenergy) :: rho
		!------------------------------------- GREEN PARAMETERS:
		COMPLEX(8), PARAMETER :: eta=DCMPLX(0.d0,0.05d0)
		COMPLEX(8) :: g !função de Green
		!------------------------------------- GAUSSIAN PARAMETERS:
		REAL(8) :: alpha1dw
		REAL(8) :: alpha2dw
		REAL(8) :: alpha1up
		REAL(8) :: alpha2up
		REAL(8) :: alpha_sum
		REAL(8), PARAMETER :: smearing=0.1d0
		REAL(8) :: factor
		REAL(8) :: arg_rho
		
		!------------------------------------- MU PARAMETERS:
		REAL(8) :: mu
		REAL(8) :: ef
		REAL(8), PARAMETER :: muconv =1.d-6
		REAL(8) :: occ_var
		REAL(8) :: integral_ne
		
		!------------------------------------- OCCUPATION PARAMETERS:
		REAL(8), PARAMETER :: CONocc=1.d-4 ! stop self-consistent occopation
		REAL(8) :: n1up
		REAL(8) :: n1dw
		REAL(8) :: n2up
		REAL(8) :: n2dw
		REAL(8) :: NEW_n1up=0.d0
		REAL(8) :: NEW_n1dw=0.d0
		REAL(8) :: NEW_n2up=0.d0
		REAL(8) :: NEW_n2dw=0.d0
		REAL(8) :: f1dw, f2dw, f1up, f2up, fsum
		INTEGER(8), PARAMETER :: DIMPHASE=3
		INTEGER(8), DIMENSION(DIMPHASE) :: phase
		INTEGER(8) :: iphase
		REAL(8), DIMENSION(1:DIMPHASE) :: convn1up
		REAL(8), DIMENSION(1:DIMPHASE) :: convn1dw
		REAL(8), DIMENSION(1:DIMPHASE) :: convn2up
		REAL(8), DIMENSION(1:DIMPHASE) :: convn2dw
		!------------------------------------- ENTROPY PARAMETERS:
		
		REAL(8), DIMENSION(1:DIMPHASE) :: TS !tempeture*entropy
		REAL(8) :: stemp !VARIABLE AUXILIAR
		REAL(8) :: s1dw !VARIABLE AUXILIAR
		REAL(8) :: s2dw !VARIABLE AUXILIAR
		REAL(8) :: s1up !VARIABLE AUXILIAR
		REAL(8) :: s2up !VARIABLE AUXILIAR
		REAL(8), DIMENSION(1:DIMPHASE) :: ET !total energy
		REAL(8) :: etemp !VARIABLE AUXILIAR
		REAL(8), DIMENSION(1:DIMPHASE) :: FREE_ENERGY !Free-energy
		
		!------------------------------------- PHASE CONVERGENCE PARAMETERS:
		!REAL(8), DIMENSION(0,NVALU*NVALne) :: TRAN_FERRO_ANTI_U
		!REAL(8), DIMENSION(0,NVALU*NVALne) :: TRAN_ANTI_PARA_U
		!REAL(8), DIMENSION(0,NVALU*NVALne) :: TRAN_FERRO_PARA_U
		!REAL(8), DIMENSION(0,NVALU*NVALne) :: TRAN_FERRO_ANTI_ne
		!REAL(8), DIMENSION(0,NVALU*NVALne) :: TRAN_ANTI_PARA_ne
		!REAL(8), DIMENSION(0,NVALU*NVALne) :: TRAN_FERRO_PARA_ne
		
		REAL(8) :: PHASE_FERRO_U !arquivo com os apenas valores de U - ferro
		REAL(8) :: PHASE_FERRO_ne !arquivo com os apenas valores de ne - ferro
		REAL(8) :: PHASE_ANTI_U !arquivo com os apenas valores de U - anti
		REAL(8) :: PHASE_ANTI_ne !arquivo com os apenas valores de ne - anti
		REAL(8) :: PHASE_PARA_U !arquivo com os apenas valores de U - para
		REAL(8) :: PHASE_PARA_ne !arquivo com os apenas valores de ne - para
		
		REAL(8), DIMENSION(edgeinix:edgefimx,edgeiniy:edgefimy,DIMHILB,DIMPHASE) :: edw_save
		REAL(8), DIMENSION(edgeinix:edgefimx,edgeiniy:edgefimy,DIMHILB,DIMPHASE) :: eup_save

END MODULE interface
