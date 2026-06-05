	real*8 function peegamma(vertex,main)

* Purpose:
* This function determines the kinematics in the PHOTON-NUCLEON center of mass
* frame, calculates some of the kinematical variables (s,t, and CM quantities
* in the 'main' structure), and returns the pion cross section.
*
*   output:
*	peegamma		!d5sigma/dEe'dOmegae'Omega	(microbarn/MeV/sr^2)

	USE structureModule
	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass, not the photon-NUCLEUS!  The model gives
* the cross section in the photon-nucleon center of mass frame.

	real*8 Q2gev,xb,thpqdeg,phpqdeg,s
	real*8 sigma_eegamma,dvcs_xsec
	real*8 FF,W1,W2,Wnuc,Wp,qmu4mp,w1p,w2p,FFrat,GEP,GMP,GEneut,GMN,w1d,w2d,tau
	real*8 tgev,FDEL
	

	! Variables calculated in transformation to gamma-NUCLEON center of mass.
        real*8 gstar,bstar,bstarx,bstary,bstarz		!beta of boost to C.M.
        real*8 nustar,qstar,qstarx,qstary,qstarz	!q in C.M.
        real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz		!p_hadron in C.M.
        real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
        real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz	!p_fermi in C.M.
        real*8 thetacm,phicm,phiqn,jacobian,jac_old,jacobian_dis

	call transform_to_cm(vertex,main,
     &		gstar,bstar,bstarx,bstary,bstarz,
     &		nustar,qstar,qstarx,qstary,qstarz,
     &		epicm,ppicm,ppicmx,ppicmy,ppicmz,
     &		ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz,
     &		etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz,
     &		thetacm,phicm,phiqn,jacobian,jac_old)

	main%thetacm = thetacm
	main%phicm = phicm
	main%pcm = ppicm
	main%davejac = jacobian
	main%johnjac = jac_old	!approx. assuming collinear boost.
c	write (6,*) jacobian,jac_old,100.*(jacobian-jac_old)/jacobian,'%'
	
	Q2gev=vertex%q2/1.0d6
	xb=vertex%q2/2.0/Mp/vertex%nu
	thpqdeg=main%theta_pq*degrad
	phpqdeg=main%phi_pq*degrad

	sigma_eegamma=dvcs_xsec(Q2gev,xb,thpqdeg,phpqdeg)/1.0d15	! dsigma/(dQ2 dxb dt dphi) in nb/GeV4--> ub/MeV4

	ntup%sigcm=sigma_eegamma*1.0d15

	s=vertex%Q2*vertex%Ein/(xb*vertex%nu)+Mp**2+Me**2
	jacobian_dis = (s-Mp**2)*xb*vertex%e%E/(2.0*pi*Mp*vertex%nu)

	tgev=main%t/1.0E6
	FFrat=1.0
	if(which_dvcs.ge.10) then
           call NFORM_XEM(14,tgev,GEP,GEneut,GMP,GMN)
c	   call fofa_best_fit(-main%t/hbarc**2,GE,GM)

	   FF  = FDEL(tgev)
	   TAU = main%t/4./Mp2  
           W1d  = FF**2*TAU*.6667*(GMN+GMP)**2                           
           W2d  = W1d+(FF*(1.0*(GEneut+TAU*GMN)+GEP+TAU*GMP)/(1.+TAU))**2   
	   Wnuc = W2d+2.0*W1d*(tan(vertex%e%theta))**2

	   qmu4mp = main%t/4./Mp2
	   W1p = GMP**2*qmu4mp
	   W2p = (GEP**2+GMP**2*qmu4mp)/(1.0+qmu4mp)
	   Wp = W2p + 2.*W1p*tan(vertex%e%theta/2.)**2
	   FFrat = Wnuc/Wp
	endif
	peegamma=FFrat*sigma_eegamma*jacobian*jacobian_dis
c	write(6,*) 'cheesy poofs',main%t/1.0d6,Wnuc,Wp,FFrat

	return
	end
	

	real*8 function dvcs_xsec(Q2,xb,thpq,phipq)
c returns cross section, dsigma/(dQ2 dxb dt dphi) in nb/GeV2
	
	implicit none
	
	real*8 Q2,xb,t,thpq,phipq
	real*8 Q2tab(46), Q2hi,Q2lo,dQ2,Q2d
	real*8 xtab(56), xhi,xlo,dx,xd
	real*8 phitab(36), phihi,philo,dphi,phid,phimax
	real*8 thtab(21), thhi,thlo,dth,thd
	real*8 c000,c001,c010,c100,c011,c101,c110,c111
	real*8 c00,c01,c10,c11
	real*8 c0,c1,c
	real*8 q2dum,xdum,tdum,phidum,thdum

	real*8 xsec(46,56,21,37)

	integer i,j,k,l
	logical first
	data first/.true./
	character*80 filename,line
	
	save

c	write(6,*) 'top of dvcs_xsec',Q2,xb,thpq,phipq
	if (first) then
	   first=.false.
	   filename='xs_table_10.6_HMS_GeV_corr.csv'
	   write(6,*) 'Initializing and reading in DVCS cross section table',filename
	   do i=1,46
	      Q2tab(i)=1.0+(i-1)*0.2
	   enddo
	
	   do i=1,56
	      xtab(i)=0.15+(i-1)*0.01
	   enddo

	   do i=1,21
	      thtab(i) = (i-1)*0.5
	   enddo

	   do i=1,37
	      phitab(i) = (i-1)*10.0
	   enddo

	   open(unit=10, file=filename, status='old')
	   read(10,'(a)') line !skip first line
	   
	   do i=1,46 ! Q2
	     do j=1,56		!xBj
		do k=1,21	!thpq
		   do l=1,36	!phipq
		      read(10,*) q2dum,xdum,tdum,phidum,thdum,xsec(i,j,k,l)
		      if(l.eq.36) then
			 xsec(i,j,k,37) = xsec(i,j,k,1)
c                        write(6,*) q2dum,xdum,tdum,phidum,thdum,xsec(i,j,k,l)
		      endif
		   enddo
		enddo
	      enddo
	   enddo
	endif

C loop over variables, find table entries
	c=0
	do i=1,46		!Q2
	   if(Q2.gt.Q2tab(i) .and. Q2.le.Q2tab(i+1)) then
	      Q2lo=Q2tab(i)
	      Q2hi=Q2tab(i+1)
	      dQ2=Q2hi-Q2lo
	      do j=1,56		!xBj
		 if(xb.gt.xtab(j) .and. xb.le.xtab(j+1)) then
		    xlo=xtab(j)
		    xhi=xtab(j+1)
		    dx=xhi-xlo
		    do k=1,21	!thpq
		       if(thpq.gt.thtab(k) .and. thpq.le.thtab(k+1)) then
			  thlo=thtab(k)
			  thhi=thtab(k+1)
			  dth=thhi-thlo
			  do l=1,37 !phipq
			     if(phipq.gt.phitab(l) .and. phipq.le.phitab(l+1)) then
				philo=phitab(l)
				phihi=phitab(l+1)
				dphi=phihi-philo

				Q2d=(Q2-Q2lo)/dQ2
				xd=(xb-xlo)/dx
				phid=(phipq-philo)/dphi
				thd=(thpq-thlo)/dth

				c000=xsec(i,j,k,l)*(1-Q2d)+xsec(i+1,j,k,l)*Q2d
				c001=xsec(i,j,k,l+1)*(1-Q2d)+xsec(i+1,j,k,l+1)*Q2d
				c010=xsec(i,j,k+1,l)*(1-Q2d)+xsec(i+1,j,k+1,l)*Q2d
				c100=xsec(i,j+1,k,l)*(1-Q2d)+xsec(i+1,j+1,k,l)*Q2d

				c011=xsec(i,j,k+1,l+1)*(1-Q2d)+xsec(i+1,j,k+1,l+1)*Q2d
				c101=xsec(i,j+1,k,l+1)*(1-Q2d)+xsec(i+1,j+1,k,l+1)*Q2d
				c110=xsec(i,j+1,k+1,l)*(1-Q2d)+xsec(i+1,j+1,k+1,l)*Q2d
				c111=xsec(i,j+1,k+1,l+1)*(1-Q2d)+xsec(i+1,j+1,k+1,l+1)*Q2d

				c00 = c000*(1-xd)+c100*xd !(0000 and 1000) and (0100 and 1100)
				c01 = c001*(1-xd)+c101*xd !(0001 and 1001) and (0101 and 1101)
				c10 = c010*(1-xd)+c110*xd !(0010 and 1010) and (0110 and 1110)
				c11 = c011*(1-xd)+c111*xd !(0011 and 1011) and (0111 and 1111)

				c0 = c00*(1-phid) + c10*phid
				c1 = c01*(1-phid) + c11*phid

				c = c0*(1-thd) + c1*thd
			     endif ! phi
			  enddo ! phi
		       endif	! theta
		    enddo	! theta
		 endif		!xb
	      enddo		!xb
	   endif		! Q2
	enddo			! Q2

	dvcs_xsec=c

	return
	end
