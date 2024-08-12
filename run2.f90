PROGRAM maestro

	USE INTERFACE

	IMPLICIT NONE
	
	!------------------------------------- DIAGONALIZATION PARAMETERS:
	INTEGER(4)            :: INFOdw
	INTEGER(4)            :: INFOup
	REAL(8), ALLOCATABLE  :: RWORKdw(:)
	REAL(8), ALLOCATABLE  :: RWORKup(:)
	COMPLEX(8),ALLOCATABLE:: WORKdw(:)
	COMPLEX(8),ALLOCATABLE:: WORKup(:)
	INTEGER(8),ALLOCATABLE:: LWORKdw(:)
	INTEGER(8),ALLOCATABLE:: LWORKup(:)
	REAL(8), ALLOCATABLE  :: e_temp(:)
	REAL(8), DIMENSION(0:Nenergy) :: integral
	integer, parameter :: ninte=20
	integer :: p
	INTEGER :: pf=1
	INTEGER :: pa=1
	INTEGER :: pp=1
	REAL(8), DIMENSION(0:Nx,0:Ny,DIMHILB) :: var
	
	INTEGER(8) :: inx !indice do vetor k em x
	INTEGER(8) :: iny !indice do vetor k em y
	real(8) :: unum

	!------------------------------------- 	INTERFACE:
	PHASE_FERRO_U=-1
	PHASE_FERRO_ne=-1
	PHASE_ANTI_U=-1
	PHASE_ANTI_ne=-1
	PHASE_PARA_U=-1
	PHASE_PARA_ne=-1
	PRINT*, kb*temperature
	print*, Nk
	BETA = 1.D0/(kb*temperature)

	!-------------------------------------  PREPARING STEPS OF ne:
	stepne = (NEMAX-ne0)/(NVALne)
	!-------------------------------------  PREPARING STEPS OF tU:
	steptU=(tUMAX-tU0)/(NVALtU)
	!-------------------------------------  DEFINE PHASE:
	phase(1)=1 !caso ferro
	phase(2)=2 !caso antiferro
	phase(3)=3 !caso para
	!DO itu=0,NVALtU
		tU=0.2d0!tU0+itu*steptU
		U=t/tU
		!DO ine=0,NVALne
			ne=1.6d0!ne0+ine*stepne
			delta= 0.05d0!0.1*ne/4
			PRINT*, "-------------------"
			PRINT*, "START PROGRAM/LOOP WITH:"
			PRINT*, "U = ", U
			PRINT*, "ne = ", ne
			PRINT*, "t=", t
			PRINT*, "temperature=", temperature
			PRINT*, "delta=", delta
			FREE_ENERGY=zero !Garantindo que Free energy inicie em zero
			ET=zero!Garantindo que total energy inicie em zero
			TS=zero!Garantindo que entropy inicie em zero
			DO iphase=1,DIMPHASE
				initial=1
				IF (phase(iphase)==1) THEN !FERRO
					NEW_n1dw=zero
					NEW_n2dw=zero
					NEW_n1up=zero
					NEW_n2up=zero
					n1dw=ne/4.d0-delta
					n2dw=ne/4.d0-delta
					n1up=ne/4.d0+delta
					n2up=ne/4.d0+delta
				ELSE IF (phase(iphase)==2) THEN !ANTI
					NEW_n1dw=zero
					NEW_n2dw=zero
					NEW_n1up=zero
					NEW_n2up=zero
					n1dw=ne/4.d0-delta
					n2dw=ne/4.d0+delta
					n1up=ne/4.d0+delta
					n2up=ne/4.d0-delta
				ELSE IF (phase(iphase)==3) THEN !PARA
					NEW_n1dw=zero
					NEW_n2dw=zero
					NEW_n1up=zero
					NEW_n2up=zero
					n1dw=ne/4.d0
					n2dw=ne/4.d0
					n1up=ne/4.d0
					n2up=ne/4.d0
				ENDIF
				PRINT*, "SEED PHASE =", phase(iphase), "1-FERRO, 2-ANTI, 3-PARA"
				PRINT*, "n1dw = ", n1dw
				PRINT*, "n2dw = ", n2dw
				PRINT*, "n1up = ", n1up
				PRINT*, "n2up = ", n2up
				PRINT*, "-------------------"
				DO WHILE (ABS(n1dw-NEW_n1dw)>CONocc.OR.ABS(n2dw-NEW_n2dw)>CONocc.OR.ABS(n1up-NEW_n1up)>CONocc.OR.ABS(n2up-NEW_n2up)>CONocc)
					PRINT*, "///////////////////////////////////////"
					PRINT*, "START WHILE FOR CONVERGENCE OCCUPATION:"
					IF (initial==1) THEN
						n1dw=n1dw
						n2dw=n2dw
						n1up=n1up
						n2up=n2up
					ELSE
						n1dw=NEW_n1dw
						n2dw=NEW_n2dw
						n1up=NEW_n1up
						n2up=NEW_n2up
					END IF
						PRINT*, "START DIAGONALIZATION..."
						DO iny= edgeiniy,edgefimy
							DO inx= edgeinix, edgefimx
								veck(:)= (/ (inx*pi)/Nx,(iny*pi)/Ny /)
								gama=-1-EXP(-i* dot_product(veck,b1+b2)) - EXP(-i* dot_product(veck,b1))- EXP(-i* dot_product(veck,b2))
								hdw(1,1)=U*n1up-U*(n1up*n1dw+n2up*n2dw)
								hdw(1,2)=t*gama
								hdw(2,1)=t*DCONJG(gama)
								hdw(2,2)=U*n2up-U*(n1up*n1dw+n2up*n2dw)
								
								hup(1,1)=U*n1dw-U*(n1up*n1dw+n2up*n2dw)
								hup(1,2)=t*gama
								hup(2,1)=t*DCONJG(gama)
								hup(2,2)=U*n2dw-U*(n1up*n1dw+n2up*n2dw)
								
								ALLOCATE(WORKdw(2*DIMHILB-1))
								ALLOCATE(LWORKdw(2*DIMHILB-1))
								ALLOCATE(RWORKdw(3*DIMHILB-2))
								ALLOCATE(e_temp(DIMHILB))
								
								CALL zheev('V','L',DIMHILB,hdw,DIMHILB,e_temp,WORKdw,2*DIMHILB-1,RWORKdw,INFOdw)
								edw(inx,iny,1)=e_temp(1)
								edw(inx,iny,2)=e_temp(2)
								
								!edw_save(inx,iny,1,iphase)=e_temp(1)
								!edw_save(inx,iny,2,iphase)=e_temp(2)
								
								vdw(inx,iny,1,1)=hdw(1,1)
								vdw(inx,iny,1,2)=hdw(2,1)
								vdw(inx,iny,2,1)=hdw(1,2)
								vdw(inx,iny,2,2)=hdw(2,2)
								
								DEALLOCATE(WORKdw)
								DEALLOCATE(LWORKdw)
								DEALLOCATE(RWORKdw)
								DEALLOCATE(e_temp)
								
								ALLOCATE(WORKup(2*DIMHILB-1))
								ALLOCATE(LWORKup(2*DIMHILB-1))
								ALLOCATE(RWORKup(3*DIMHILB-2))
								ALLOCATE(e_temp(DIMHILB))
								
								CALL zheev('V','L',DIMHILB,hup,DIMHILB,e_temp,WORKup,2*DIMHILB-1,RWORKup,INFOup)
								eup(inx,iny,1)=e_temp(1)
								eup(inx,iny,2)=e_temp(2)
								
								!eup_save(inx,iny,1,iphase)=e_temp(1)
								!eup_save(inx,iny,2,iphase)=e_temp(2)
								
								vup(inx,iny,1,1)=hup(1,1)
								vup(inx,iny,1,2)=hup(2,1)
								vup(inx,iny,2,1)=hup(1,2)
								vup(inx,iny,2,2)=hup(2,2)

								DEALLOCATE(WORKup)
								DEALLOCATE(LWORKup)
								DEALLOCATE(RWORKup)
								DEALLOCATE(e_temp)
								
							END DO
						END DO
						PRINT*, "///////////////////////////////////////"
						! definindo os valores da energia:
						minenergy(1)=MINVAL(edw)
						minenergy(2)=MINVAL(eup)
						minenergy=MINVAL(minenergy)
						maxenergy(1)=MAXVAL(edw)
						maxenergy(2)=MAXVAL(eup)
						maxenergy=MAXVAL(maxenergy)
						step_energy=val_step_energy
						PRINT*, "step energy=", step_energy
						nsteps=INT(((maxenergy(1)+0.5d0)-(minenergy(1)-0.5d0))/val_step_energy)
						PRINT*, "NSTEPS=", nsteps
						!construindo rho: VIA GREEN'S FUNCTION:
						
						!DO ienergy=0,nsteps
						!	w(ienergy)=minenergy(1)-0.5d0+ienergy*step_energy
						!	g=DCMPLX(0.d0,0.d0)
						!	DO iny= 1,Ny
						!		DO inx= 1, Nx
						!			g=g+1/(w(ienergy)+eta-edw(inx,iny,1)) + 1/(w(ienergy)+eta-edw(inx,iny,2))+&
						!				& 1/(w(ienergy)+eta-eup(inx,iny,1)) + 1/(w(ienergy)+eta-eup(inx,iny,2))
						!		ENDDO
						!	ENDDO
						!	rho(ienergy)=-DIMAG(g)/(Nk*pi) 
						!ENDDO
						
						!construindo rho: VIA GAUSSIAN FUNCTION:
						factor=1.d0/(smearing*REAL(Nk)*DSQRT(twopi/2.d0))
						PRINT*, "factor", factor
						PRINT*, "smearing", smearing
						PRINT*, "smearing2", smearing**2.d0
						rho=zero
						DO ienergy=0,nsteps
							w(ienergy)=minenergy(1)-0.5d0+ienergy*step_energy
							arg_rho=zero
							DO iny= edgeiniy,edgeborday
								DO inx= edgeinix, edgebordax
									alpha1dw=DEXP(-1.d0*((w(ienergy)-edw(inx,iny,1))**2.d0)/(smearing**2.d0))
									alpha2dw=DEXP(-1.d0*((w(ienergy)-edw(inx,iny,2))**2.d0)/(smearing**2.d0))
									alpha1up=DEXP(-1.d0*((w(ienergy)-eup(inx,iny,1))**2.d0)/(smearing**2.d0))
									alpha2up=DEXP(-1.d0*((w(ienergy)-eup(inx,iny,2))**2.d0)/(smearing**2.d0))
									alpha_sum = alpha1dw+alpha2dw+alpha1up+alpha2up
									arg_rho= arg_rho + alpha_sum
								END DO
							END DO
							rho(ienergy)=factor*arg_rho
						END DO
						!calculando a integral de rho para encontrar ef: 
						integral = 0.d0
						DO ienergy=1,nsteps-1
							integral(ienergy) = integral(ienergy-1) + (rho(ienergy)+rho(ienergy+1))
						END DO
						integral=integral*step_energy/2.d0
						print*, "integral=",integral(nsteps-1)
						PRINT*, "============================"
						PRINT*, "MINENERGY = ", minenergy(1)
						PRINT*, "MAXENERGY = ", maxenergy(1)
						PRINT*, "START EF EVALUT ... "
						DO ienergy=0,nsteps
							IF (integral(ienergy)>=ne) THEN
								EXIT
							END IF
						END DO
						PRINT*, INTEGRAL(IENERGY)
						PRINT*, IENERGY
						!DEFINE THE EF:
						IF (integral(ienergy)==ne) THEN
							ef=w(ienergy)
						ELSE
							ef=w(ienergy-1)+step_energy*(ne-integral(ienergy-1))/                       &
								&(integral(ienergy)-integral(ienergy-1))
						END IF
						PRINT*, "EF = ", EF
						PRINT*, "============================"
						! Calculo de mu:
						call mueval
						PRINT*, "mu=", mu
						!Calculando a integral com mu convergido:
						integral_ne = 0.d0
						DO iny= edgeiniy,edgeborday
								DO inx= edgeinix, edgebordax
									occ_var =  1.D0/(DEXP(BETA*(edw(inx,iny,1)-mu)) + 1.D0)+&
											& 1.D0/(DEXP(BETA*(edw(inx,iny,2)-mu)) + 1.D0)+&
											& 1.D0/(DEXP(BETA*(eup(inx,iny,1)-mu)) + 1.D0)+&
											& 1.D0/(DEXP(BETA*(eup(inx,iny,2)-mu)) + 1.D0)
									integral_ne=integral_ne+occ_var
								END DO
						END DO
						integral_ne = integral_ne/(Nk)
						PRINT*, 'integral c/ mu convergido=', integral_ne, 'ne=', ne
						! Calculo das novas ocupações
						NEW_n1up=0.0d0
						NEW_n1dw=0.0d0
						NEW_n2up=0.0d0
						NEW_n2dw=0.0d0
						DO iny= edgeiniy,edgeborday
							DO inx= edgeinix, edgebordax
								f1dw = 1.D0/(DEXP(BETA*(edw(inx,iny,1)-mu)) + 1.D0)
								f2dw = 1.D0/(DEXP(BETA*(edw(inx,iny,2)-mu)) + 1.D0)
								f1up = 1.D0/(DEXP(BETA*(eup(inx,iny,1)-mu)) + 1.D0)
								f2up = 1.D0/(DEXP(BETA*(eup(inx,iny,2)-mu)) + 1.D0)
								
								NEW_n1dw = NEW_n1dw+f1dw*DREAL(vdw(inx,iny,1,1)*DCONJG(vdw(inx,iny,1,1)))/(Nk) &
											&+ f2dw*DREAL(vdw(inx,iny,1,2)*DCONJG(vdw(inx,iny,1,2)))/(Nk)
								NEW_n2dw = NEW_n2dw+f1dw*DREAL(vdw(inx,iny,2,1)*DCONJG(vdw(inx,iny,2,1)))/(Nk) &
											&+ f2dw*DREAL(vdw(inx,iny,2,2)*DCONJG(vdw(inx,iny,2,2)))/(Nk)
								NEW_n1up = NEW_n1up+f1up*DREAL(vup(inx,iny,1,1)*DCONJG(vup(inx,iny,1,1)))/(Nk) &
											&+ f2up*DREAL(vup(inx,iny,1,2)*DCONJG(vup(inx,iny,1,2)))/(Nk)
								NEW_n2up = NEW_n2up+f1up*DREAL(vup(inx,iny,2,1)*DCONJG(vup(inx,iny,2,1)))/(Nk) &
											&+ f2up*DREAL(vup(inx,iny,2,2)*DCONJG(vup(inx,iny,2,2)))/(Nk)
							ENDDO
						ENDDO
						initial=0
						PRINT *,"newn1dw=", NEW_n1dw
						PRINT *,"newn1up=", NEW_n1up
						PRINT *,"newn2dw=", NEW_n2dw
						PRINT *,"newn2up=", NEW_n2up
						PRINT *,"soma occ=", NEW_n1dw+NEW_n1up+NEW_n2dw+ NEW_n2up
				ENDDO !Fim do while para convergência das ocupações
				PRINT *, 'phase=', phase(iphase)
				convn1up(iphase) = NEW_n1up
				convn2up(iphase) = NEW_n2up
				convn1dw(iphase) = NEW_n1dw
				convn2dw(iphase) = NEW_n2dw
				PRINT*, "***************************"
				PRINT*, "CONVERGED PHASE:"
				PRINT *,"newn1dw=", convn1dw(iphase)
				PRINT *,"newn1up=", convn1up(iphase)
				PRINT *,"newn2dw=", convn2dw(iphase)
				PRINT *,"newn2up=", convn2up(iphase)
				PRINT*, "***************************"
				!DIAGONALIZANDO NOVAMENTE, JÁ COM AS OCUPAÇÕES CONVERGIDAS:
				PRINT*, "STARTING DIAGONALIZATION WITH OCC CONVERGED CALCULATION"
				unum=-U*(convn1dw(iphase)*convn1up(iphase)+convn2dw(iphase)*convn2up(iphase))
				edw=zero
				eup=zero
				DO iny= edgeiniy,edgefimy
					DO inx= edgeinix, edgefimx
						veck(:)= (/ (inx*pi)/Nx,(iny*pi)/Ny /)
						gama=-1-EXP(-i* dot_product(veck,b1+b2)) - EXP(-i* dot_product(veck,b1))- EXP(-i* dot_product(veck,b2))
						hdw(1,1)=U*convn1up(iphase)+unum
						hdw(1,2)=t*gama
						hdw(2,1)=t*DCONJG(gama)
						hdw(2,2)=U*convn2up(iphase)+unum
						
						hup(1,1)=U*convn1dw(iphase)+unum
						hup(1,2)=t*gama
						hup(2,1)=t*DCONJG(gama)
						hup(2,2)=U*convn2dw(iphase)+unum
								
						ALLOCATE(WORKdw(2*DIMHILB-1))
						ALLOCATE(LWORKdw(2*DIMHILB-1))
						ALLOCATE(RWORKdw(3*DIMHILB-2))
						ALLOCATE(e_temp(DIMHILB))
								
						CALL zheev('V','L',DIMHILB,hdw,DIMHILB,e_temp,WORKdw,2*DIMHILB-1,RWORKdw,INFOdw)
						edw(inx,iny,1)=e_temp(1)
						edw(inx,iny,2)=e_temp(2)

						vdw(inx,iny,1,1)=hdw(1,1)
						vdw(inx,iny,1,2)=hdw(2,1)
						vdw(inx,iny,2,1)=hdw(1,2)
						vdw(inx,iny,2,2)=hdw(2,2)
								
						DEALLOCATE(WORKdw)
						DEALLOCATE(LWORKdw)
						DEALLOCATE(RWORKdw)
						DEALLOCATE(e_temp)
								
						ALLOCATE(WORKup(2*DIMHILB-1))
						ALLOCATE(LWORKup(2*DIMHILB-1))
						ALLOCATE(RWORKup(3*DIMHILB-2))
						ALLOCATE(e_temp(DIMHILB))
								
						CALL zheev('V','L',DIMHILB,hup,DIMHILB,e_temp,WORKup,2*DIMHILB-1,RWORKup,INFOup)
						eup(inx,iny,1)=e_temp(1)
						eup(inx,iny,2)=e_temp(2)

						vup(inx,iny,1,1)=hup(1,1)
						vup(inx,iny,1,2)=hup(2,1)
						vup(inx,iny,2,1)=hup(1,2)
						vup(inx,iny,2,2)=hup(2,2)

						DEALLOCATE(WORKup)
						DEALLOCATE(LWORKup)
						DEALLOCATE(RWORKup)
						DEALLOCATE(e_temp)
								
					END DO
				END DO
				
				minenergy(1)=MINVAL(edw)
				minenergy(2)=MINVAL(eup)
				minenergy=MINVAL(minenergy)
				maxenergy(1)=MAXVAL(edw)
				maxenergy(2)=MAXVAL(eup)
				maxenergy=MAXVAL(maxenergy)
				step_energy=val_step_energy
				nsteps=INT(((maxenergy(1)+0.5d0)-(minenergy(1)-0.5d0))/val_step_energy)
				
				factor=1.d0/(smearing*REAL(Nk)*DSQRT(twopi/2.d0))
				rho=zero
				DO ienergy=0,nsteps
					w(ienergy)=minenergy(1)-0.5d0+ienergy*step_energy
					arg_rho=zero
					DO iny= edgeiniy,edgeborday
						DO inx= edgeinix, edgebordax
							alpha1dw=DEXP(-1.d0*((w(ienergy)-edw(inx,iny,1))**2.d0)/(smearing**2.d0))
							alpha2dw=DEXP(-1.d0*((w(ienergy)-edw(inx,iny,2))**2.d0)/(smearing**2.d0))
							alpha1up=DEXP(-1.d0*((w(ienergy)-eup(inx,iny,1))**2.d0)/(smearing**2.d0))
							alpha2up=DEXP(-1.d0*((w(ienergy)-eup(inx,iny,2))**2.d0)/(smearing**2.d0))
							alpha_sum = alpha1dw+alpha2dw+alpha1up+alpha2up
							arg_rho= arg_rho + alpha_sum
						END DO
					END DO
					rho(ienergy)=factor*arg_rho
				END DO
				!calculando a integral de rho para encontrar ef: 
				integral = 0.d0
				DO ienergy=1,nsteps-1
					integral(ienergy) = integral(ienergy-1) + (rho(ienergy)+rho(ienergy+1))
				END DO
				integral=integral*step_energy/2.d0
				DO ienergy=0,nsteps
					IF (integral(ienergy)>=ne) THEN
						EXIT
					END IF
				END DO
				!DEFINE THE EF:
				IF (integral(ienergy)==ne) THEN
					ef=w(ienergy)
				ELSE
					ef=w(ienergy-1)+step_energy*(ne-integral(ienergy-1))/                       &
						&(integral(ienergy)-integral(ienergy-1))
				END IF
				! Calculo de mu:
				call mueval
				PRINT*, "mu=", mu
				PRINT*, "STARTING TOTAL ENERGY AND ENTROPY CALCULATION"
				stemp=0.0d0
				s1dw=0.d0
				s2dw=0.d0
				s1up=0.d0
				s2up=0.d0
				etemp=0.0d0
				DO iny= edgeiniy,edgeborday
					DO inx= edgeinix, edgebordax
						f1dw = 1.D0/(DEXP(BETA*(edw(inx,iny,1)-mu)) + 1.D0)
						f2dw = 1.D0/(DEXP(BETA*(edw(inx,iny,2)-mu)) + 1.D0)
						f1up = 1.D0/(DEXP(BETA*(eup(inx,iny,1)-mu)) + 1.D0)
						f2up = 1.D0/(DEXP(BETA*(eup(inx,iny,2)-mu)) + 1.D0)
						
						etemp = etemp+ f1dw*edw(inx,iny,1) + f2dw*edw(inx,iny,2)+&
									 & f1up*eup(inx,iny,1) + f2up*eup(inx,iny,2)
						
						IF ((f1dw>1.d-100).AND.(1.d0-f1dw>1.d-100)) THEN
							s1dw=s1dw+(DLOG(f1dw)*f1dw+ &
							&DLOG(1.d0-f1dw)*(1.d0-f1dw))
						END IF	
						IF ((f2dw>1.d-100).AND.(1.d0-f2dw>1.d-100)) THEN
							s2dw=s2dw+(DLOG(f2dw)*f2dw+ &
							&DLOG(1.d0-f2dw)*(1.d0-f2dw))
						END IF	
						IF ((f1up>1.d-100).AND.(1.d0-f1up>1.d-100)) THEN
							s1up=s1up+(DLOG(f1up)*f1up+ &
							&DLOG(1.d0-f1up)*(1.d0-f1up))
						END IF	
						IF ((f2up>1.d-100).AND.(1.d0-f2up>1.d-100)) THEN
							s2up=s2up+(DLOG(f2up)*f2up+ &
							&DLOG(1.d0-f2up)*(1.d0-f2up))
						END IF		  
					END DO
				END DO
				ET(iphase) = etemp/(Nk)
				stemp=s1dw+s2dw+s1up+s2up
				TS(iphase) = -stemp/(Nk*BETA) 
				FREE_ENERGY(iphase) = ET(iphase) - TS(iphase)
				PRINT*, 'TOTAL ENERGY = ', ET
				PRINT*, 'TEMPERATURE*ENTROPY = ', TS
				PRINT*, 'FREE ENERGY = ', FREE_ENERGY
			END DO !Fim loop das phases
			
			!OPEN(UNIT=1,FILE='FREE_ENERGY.txt')
			!	DO iphase=1,DIMPHASE
			!		WRITE(1,*) FREE_ENERGY(iphase)
			!	END DO
			!CLOSE
			
			!OPEN(UNIT=1,FILE='ulist.txt')
			!	DO iphase=1,DIMPHASE
			!		WRITE(1,*) U
			!	END DO
			!CLOSE
			!if (abs(FREE_ENERGY(1) - FREE_ENERGY(2))<1.d-6) THEN
			!        TRAN_FERRO_ANTI_U(1)=t/U
			!        TRAN_FERRO_ANTI_ne(1)=ne/2.d0
			!END IF
			!if (abs(FREE_ENERGY(2) - FREE_ENERGY(3))<1.d-6) THEN
			!        TRAN_ANTI_PARA_U(1)=t/U
			!        TRAN_ANTI_PARA_ne(1)=ne/2.d0
			!END IF
			!if (abs(FREE_ENERGY(1) - FREE_ENERGY(3))<1.d-8) THEN
			!        TRAN_FERRO_PARA_U(1)=t/U
			!        TRAN_FERRO_PARA_ne(1)=ne/2
			!END IF
			IF (MINVAL(FREE_ENERGY)==FREE_ENERGY(1)) THEN
					PHASE_FERRO_U = U
					PHASE_FERRO_ne= ne/2
					!!PATH IN BRILLOUIN:
					!DO inx=0,Nx-1
					!	vecpath(inx,1)=(inx*pi)/Nx
					!	vecpath(inx,2)=(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+Nx,1)=pi-(inx*pi)/Nx
					!	vecpath(inx+Nx,2)=pi-(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+2*Nx,1)=0.d0 
					!	vecpath(inx+2*Nx,2)=(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+3*Nx,1)=0.d0
					!	vecpath(inx+3*Nx,2)=pi-(inx*pi)/Nx
					!END DO
					
					!n1up=convn1up(1)
					!n2up=convn2up(1)
					!n1dw=convn1dw(1)
					!n2dw=convn2dw(1)
					
					!DO inx= 0,4*Nx
					!	veck(:)= (/ vecpath(inx,1),vecpath(inx,2) /)
					!	gama=-1-EXP(-i* dot_product(veck,b1+b2)) - EXP(-i* dot_product(veck,b1))- EXP(-i* dot_product(veck,b2))
					!	hdw(1,1)=U*n1up
					!	hdw(1,2)=t*gama
					!	hdw(2,1)=t*DCONJG(gama)
					!	hdw(2,2)=U*n2up
							
					!	hup(1,1)=U*n1dw
					!	hup(1,2)=t*gama
					!	hup(2,1)=t*DCONJG(gama)
					!	hup(2,2)=U*n2dw
							
					!	ALLOCATE(WORKdw(2*DIMHILB-1))
					!	ALLOCATE(LWORKdw(2*DIMHILB-1))
					!	ALLOCATE(RWORKdw(3*DIMHILB-2))
					!	ALLOCATE(e_temp(DIMHILB))
							
					!	CALL zheev('V','L',DIMHILB,hdw,DIMHILB,e_temp,WORKdw,2*DIMHILB-1,RWORKdw,INFOdw)
					!	edw_path(inx,1)=e_temp(1)
					!	edw_path(inx,2)=e_temp(2)
								
					!	DEALLOCATE(WORKdw)
					!	DEALLOCATE(LWORKdw)
					!	DEALLOCATE(RWORKdw)
					!	DEALLOCATE(e_temp)
								
					!	ALLOCATE(WORKup(2*DIMHILB-1))
					!	ALLOCATE(LWORKup(2*DIMHILB-1))
					!	ALLOCATE(RWORKup(3*DIMHILB-2))
					!	ALLOCATE(e_temp(DIMHILB))
								
					!	CALL zheev('V','L',DIMHILB,hup,DIMHILB,e_temp,WORKup,2*DIMHILB-1,RWORKup,INFOup)
					!	eup_path(inx,1)=e_temp(1)
					!	eup_path(inx,2)=e_temp(2)

					!	DEALLOCATE(WORKup)
					!	DEALLOCATE(LWORKup)
					!	DEALLOCATE(RWORKup)
					!	DEALLOCATE(e_temp)
					!END DO
					
					!minenergy(1)=MINVAL(edw_path)
					!minenergy(2)=MINVAL(eup_path)
					!minenergy=MINVAL(minenergy)
					!maxenergy(1)=MAXVAL(edw_path)
					!maxenergy(2)=MAXVAL(eup_path)
					!maxenergy=MAXVAL(maxenergy)
					!step_energy=(maxenergy(1)+0.5-(minenergy(1)-0.5))/Nenergy
					
					!DO ienergy=0,Nenergy
					!	w(ienergy)=minenergy(1)-0.5+ienergy*step_energy
					!	g=DCMPLX(0.d0,0.d0)
					!	DO iny= 0,Ny-1
					!		DO inx= 0, Nx-1
					!			g=g+1/(w(ienergy)+eta-edw_save(inx,iny,1,1)) + 1/(w(ienergy)+eta-edw_save(inx,iny,2,1))+&
					!				& 1/(w(ienergy)+eta-eup_save(inx,iny,1,1)) + 1/(w(ienergy)+eta-eup_save(inx,iny,2,1))
					!		ENDDO
					!	ENDDO
					!	rho(ienergy)=-DIMAG(g)/(Nk*pi) 
					!ENDDO
					
					!OPEN(UNIT=1,FILE='energy_ferro.txt')
					!do inx=0,4*Nx
					!	WRITE(1,*) edw_path(inx,1), edw_path(inx,2), eup_path(inx,1), eup_path(inx,2)
					!end do
					!CLOSE(1)
					!OPEN(UNIT=1,FILE='rho_ferro.txt')
					!do ienergy=0,Nenergy
					!	WRITE(1,*) rho(ienergy), w(ienergy)
					!end do
					!CLOSE(1)
					PRINT*, "FERROMAGNETIC PHASE"
					OPEN(UNIT=1,FILE='PHASE_FERRO_U.txt',position="append")
							WRITE(1,*) PHASE_FERRO_U
					CLOSE(1)
					OPEN(UNIT=1,FILE='PHASE_FERRO_ne.txt',position="append")
							WRITE(1,*) PHASE_FERRO_ne
					CLOSE(1)

			END IF
			IF (MINVAL(FREE_ENERGY)==FREE_ENERGY(2)) THEN
					PHASE_ANTI_U = U
					PHASE_ANTI_ne= ne/2
					!!PATH IN BRILLOUIN:
					!DO inx=0,Nx-1
					!	vecpath(inx,1)=(inx*pi)/Nx
					!	vecpath(inx,2)=(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+Nx,1)=pi-(inx*pi)/Nx
					!	vecpath(inx+Nx,2)=pi+(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+2*Nx,1)=0.d0 
					!	vecpath(inx+2*Nx,2)=2*pi-(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+3*Nx,1)=0.d0
					!	vecpath(inx+3*Nx,2)=pi-(inx*pi)/Nx
					!END DO
					
					!n1up=convn1up(2)
					!n2up=convn2up(2)
					!n1dw=convn1dw(2)
					!n2dw=convn2dw(2)
					
					!DO inx= 0,4*Nx
					!	veck(:)= vecpath(inx,:)
					!	gama=-1-EXP(-i* dot_product(veck,b1+b2)) - EXP(-i* dot_product(veck,b1))- EXP(-i* dot_product(veck,b2))
					!	hdw(1,1)=U*n1up
					!	hdw(1,2)=t*gama
					!	hdw(2,1)=t*DCONJG(gama)
					!	hdw(2,2)=U*n2up
							
					!	hup(1,1)=U*n1dw
					!	hup(1,2)=t*gama
					!	hup(2,1)=t*DCONJG(gama)
					!	hup(2,2)=U*n2dw
							
					!	ALLOCATE(WORKdw(2*DIMHILB-1))
					!	ALLOCATE(LWORKdw(2*DIMHILB-1))
					!	ALLOCATE(RWORKdw(3*DIMHILB-2))
					!	ALLOCATE(e_temp(DIMHILB))
							
					!	CALL zheev('V','L',DIMHILB,hdw,DIMHILB,e_temp,WORKdw,2*DIMHILB-1,RWORKdw,INFOdw)
					!	edw_path(inx,1)=e_temp(1)
					!	edw_path(inx,2)=e_temp(2)
								
					!	DEALLOCATE(WORKdw)
					!	DEALLOCATE(LWORKdw)
					!	DEALLOCATE(RWORKdw)
					!	DEALLOCATE(e_temp)
								
					!	ALLOCATE(WORKup(2*DIMHILB-1))
					!	ALLOCATE(LWORKup(2*DIMHILB-1))
					!	ALLOCATE(RWORKup(3*DIMHILB-2))
					!	ALLOCATE(e_temp(DIMHILB))
								
					!	CALL zheev('V','L',DIMHILB,hup,DIMHILB,e_temp,WORKup,2*DIMHILB-1,RWORKup,INFOup)
					!	eup_path(inx,1)=e_temp(1)
					!	eup_path(inx,2)=e_temp(2)

					!	DEALLOCATE(WORKup)
					!	DEALLOCATE(LWORKup)
					!	DEALLOCATE(RWORKup)
					!	DEALLOCATE(e_temp)
					!END DO
					
					!minenergy(1)=MINVAL(edw_path)
					!minenergy(2)=MINVAL(eup_path)
					!minenergy=MINVAL(minenergy)
					!maxenergy(1)=MAXVAL(edw_path)
					!maxenergy(2)=MAXVAL(eup_path)
					!maxenergy=MAXVAL(maxenergy)
					!step_energy=(maxenergy(1)+0.5-(minenergy(1)-0.5))/Nenergy
					
					!DO ienergy=0,Nenergy
					!	w(ienergy)=minenergy(1)-0.5+ienergy*step_energy
					!	g=DCMPLX(0.d0,0.d0)
					!	DO iny= 0,Ny-1
					!		DO inx= 0, Nx-1
					!			g=g+1/(w(ienergy)+eta-edw_save(inx,iny,1,2)) + 1/(w(ienergy)+eta-edw_save(inx,iny,2,2))+&
					!				& 1/(w(ienergy)+eta-eup_save(inx,iny,1,2)) + 1/(w(ienergy)+eta-eup_save(inx,iny,2,2))
					!		ENDDO
					!	ENDDO
					!	rho(ienergy)=-DIMAG(g)/(Nk*pi) 
					!ENDDO
					
					!OPEN(UNIT=1,FILE='energy_anti.txt')
					!do inx=0,4*Nx
					!	WRITE(1,*) edw_path(inx,1), edw_path(inx,2), eup_path(inx,1), eup_path(inx,2)
					!end do
					!CLOSE(1)
					!OPEN(UNIT=1,FILE='rho_anti.txt')
					!do ienergy=0,Nenergy
					!	WRITE(1,*) rho(ienergy), w(ienergy)
					!end do
					!CLOSE(1)
					PRINT*, "ANTIFERROMAGNETIC PHASE"
					OPEN(UNIT=1,FILE='PHASE_ANTI_U.txt',position="append")
						WRITE(1,*) PHASE_ANTI_U
					CLOSE(1)
					OPEN(UNIT=1,FILE='PHASE_ANTI_ne.txt',position="append")
						WRITE(1,*) PHASE_ANTI_ne
					CLOSE(1)

			END IF
			IF (MINVAL(FREE_ENERGY)==FREE_ENERGY(3)) THEN
					PHASE_PARA_U = U
					PHASE_PARA_ne= ne/2
					!!PATH IN BRILLOUIN:
					!DO inx=0,Nx-1
					!	vecpath(inx,1)=(inx*pi)/Nx
					!	vecpath(inx,2)=(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+Nx,1)=pi-(inx*pi)/Nx
					!	vecpath(inx+Nx,2)=pi+(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+2*Nx,1)=0.d0 
					!	vecpath(inx+2*Nx,2)=2*pi-(inx*pi)/Nx
					!END DO
					!DO inx=0,Nx-1
					!	vecpath(inx+3*Nx,1)=0.d0
					!	vecpath(inx+3*Nx,2)=pi-(inx*pi)/Nx
					!END DO
					
					!n1up=convn1up(3)
					!n2up=convn2up(3)
					!n1dw=convn1dw(3)
					!n2dw=convn2dw(3)
					
					!DO inx= 0,4*Nx
					!	veck(:)= vecpath(inx,:)
					!	gama=-1-EXP(-i* dot_product(veck,b1+b2)) - EXP(-i* dot_product(veck,b1))- EXP(-i* dot_product(veck,b2))
					!	hdw(1,1)=U*n1up
					!	hdw(1,2)=t*gama
					!	hdw(2,1)=t*DCONJG(gama)
					!	hdw(2,2)=U*n2up
							
					!	hup(1,1)=U*n1dw
					!	hup(1,2)=t*gama
					!	hup(2,1)=t*DCONJG(gama)
					!	hup(2,2)=U*n2dw
							
					!	ALLOCATE(WORKdw(2*DIMHILB-1))
					!	ALLOCATE(LWORKdw(2*DIMHILB-1))
					!	ALLOCATE(RWORKdw(3*DIMHILB-2))
					!	ALLOCATE(e_temp(DIMHILB))
							
					!	CALL zheev('V','L',DIMHILB,hdw,DIMHILB,e_temp,WORKdw,2*DIMHILB-1,RWORKdw,INFOdw)
					!	edw_path(inx,1)=e_temp(1)
					!	edw_path(inx,2)=e_temp(2)
								
					!	DEALLOCATE(WORKdw)
					!	DEALLOCATE(LWORKdw)
					!	DEALLOCATE(RWORKdw)
					!	DEALLOCATE(e_temp)
								
					!	ALLOCATE(WORKup(2*DIMHILB-1))
					!	ALLOCATE(LWORKup(2*DIMHILB-1))
					!	ALLOCATE(RWORKup(3*DIMHILB-2))
					!	ALLOCATE(e_temp(DIMHILB))
								
					!	CALL zheev('V','L',DIMHILB,hup,DIMHILB,e_temp,WORKup,2*DIMHILB-1,RWORKup,INFOup)
					!	eup_path(inx,1)=e_temp(1)
					!	eup_path(inx,2)=e_temp(2)

					!	DEALLOCATE(WORKup)
					!	DEALLOCATE(LWORKup)
					!	DEALLOCATE(RWORKup)
					!	DEALLOCATE(e_temp)
					!END DO
					
					!minenergy(1)=MINVAL(edw_path)
					!minenergy(2)=MINVAL(eup_path)
					!minenergy=MINVAL(minenergy)
					!maxenergy(1)=MAXVAL(edw_path)
					!maxenergy(2)=MAXVAL(eup_path)
					!maxenergy=MAXVAL(maxenergy)
					!step_energy=(maxenergy(1)+0.5-(minenergy(1)-0.5))/Nenergy
					
					!DO ienergy=0,Nenergy
					!	w(ienergy)=minenergy(1)-0.5+ienergy*step_energy
					!	g=DCMPLX(0.d0,0.d0)
					!	DO iny= 0,Ny-1
					!		DO inx= 0, Nx-1
					!			g=g+1/(w(ienergy)+eta-edw_save(inx,iny,1,3)) + 1/(w(ienergy)+eta-edw_save(inx,iny,2,3))+&
					!				& 1/(w(ienergy)+eta-eup_save(inx,iny,1,3)) + 1/(w(ienergy)+eta-eup_save(inx,iny,2,3))
					!		ENDDO
					!	ENDDO
					!	rho(ienergy)=-DIMAG(g)/(Nk*pi) 
					!ENDDO
					
					!OPEN(UNIT=1,FILE='energy_para.txt')
					!do inx=0,4*Nx
					!	WRITE(1,*) edw_path(inx,1), edw_path(inx,2), eup_path(inx,1), eup_path(inx,2)
					!end do
					!CLOSE(1)
					!OPEN(UNIT=1,FILE='rho_para.txt')
					!do ienergy=0,Nenergy
					!	WRITE(1,*) rho(ienergy), w(ienergy)
					!end do
					!CLOSE(1)
					PRINT*, "PARAMAGNETIC PHASE"
					OPEN(UNIT=1,FILE='PHASE_PARA_U.txt',position="append")
						WRITE(1,*) PHASE_PARA_U
					CLOSE(1)
					OPEN(UNIT=1,FILE='PHASE_PARA_ne.txt',position="append")
							WRITE(1,*) PHASE_PARA_ne
					CLOSE(1)

			END IF
			
		!END DO!FIM DO PARA NE
		
	!END DO!FIM DO PARA U

END PROGRAM maestro



Subroutine mueval
	USE INTERFACE
	IMPLICIT NONE

	INTEGER :: j, q
	REAL(8) :: chp, chpminus, chpplus, f, fplus, fminus, ff
	REAL(8) :: occplus, occminus, occ
	INTEGER :: COUNTWHILEMU
	
	chpplus=ef+1.1d0
	chpminus=ef-1.1d0
	chp=ef
	fplus=0.d0
	fminus=0.d0
	DO q= edgeiniy,edgeborday
		DO j= edgeinix, edgebordax
			occplus =  1.D0/(DEXP(BETA*(edw(j,q,1)-chpplus)) + 1.D0)+&
					 & 1.D0/(DEXP(BETA*(edw(j,q,2)-chpplus)) + 1.D0)+&
					 & 1.D0/(DEXP(BETA*(eup(j,q,1)-chpplus)) + 1.D0)+&
					 & 1.D0/(DEXP(BETA*(eup(j,q,2)-chpplus)) + 1.D0)
			fplus=fplus+occplus
			occminus = 1.D0/(DEXP(BETA*(edw(j,q,1)-chpminus)) + 1.D0)+&
					 & 1.D0/(DEXP(BETA*(edw(j,q,2)-chpminus)) + 1.D0)+&
					 & 1.D0/(DEXP(BETA*(eup(j,q,1)-chpminus)) + 1.D0)+&
					 & 1.D0/(DEXP(BETA*(eup(j,q,2)-chpminus)) + 1.D0)
			fminus=fminus+occminus
		END DO
	END DO
	
	fplus=(fplus/REAL(Nk))-ne
	fminus=(fminus/REAL(Nk))-ne
	ff=fplus*fminus
	!PRINT*,"EF=", EF
	IF (ff>0) THEN
		PRINT*,'PROBLEM IN MU_EVAL: ff IS POSITIVE'
		STOP
	END IF
	COUNTWHILEMU = 0
	PRINT*, "++++++++++++++++++++++++++++"
	PRINT*, "START WHILE TO EVALUE MU ..."
	DO WHILE (DABS(chpplus-chpminus)>muconv)
		!PRINT*, "COUNT = ", COUNTWHILEMU
		f=0.d0
		fminus=0.d0
		chp=0.5d0*(chpplus+chpminus)

		DO q=edgeiniy, edgeborday
            DO j=edgeinix, edgebordax
				!Ocupação com muminus:
				occminus = 1.D0/(DEXP(BETA*(edw(j,q,1)-chpminus)) + 1.D0)+&
						 & 1.D0/(DEXP(BETA*(edw(j,q,2)-chpminus)) + 1.D0)+&
						 & 1.D0/(DEXP(BETA*(eup(j,q,1)-chpminus)) + 1.D0)+&
						 & 1.D0/(DEXP(BETA*(eup(j,q,2)-chpminus)) + 1.D0)
				fminus=fminus+occminus
				!Ocupação com mu:
                occ=  1.D0/(DEXP(BETA*(edw(j,q,1)-chp)) + 1.D0)+&
					& 1.D0/(DEXP(BETA*(edw(j,q,2)-chp)) + 1.D0)+&
					& 1.D0/(DEXP(BETA*(eup(j,q,1)-chp)) + 1.D0)+&
					& 1.D0/(DEXP(BETA*(eup(j,q,2)-chp)) + 1.D0)
                f=f+occ
            END DO   ! fim do Ky
         END DO   ! Fim do Kx

		 fminus=(fminus/REAL(Nk))-ne
		 f=(f/REAL(Nk))-ne

		 IF (fminus*f<0) THEN
			chpplus=chp
		 ELSE
			chpminus=chp
		 END IF
		 COUNTWHILEMU=COUNTWHILEMU+1
	END DO   ! Fim do while
	PRINT*, "++++++++++++++++++++++++++++"
	mu=0.5d0*(chpplus+chpminus)
end subroutine mueval
