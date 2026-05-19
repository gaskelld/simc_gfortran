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

	real*8 Q2gev,xb,thpqdeg,phpqdeg
	real*8 sigma_eegamma,dvcs_xsec

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
	main%johnjac = jac_old		!approx. assuming collinear boost.
	
	Q2gev=vertex%q2/1.0d6
	xb=vertex%q2/2.0/Mp/vertex%nu
	thpqdeg=main%theta_pq*degrad
	phpqdeg=main%phi_pq*degrad

	sigma_eegamma=dvcs_xsec(Q2gev,xb,thpqdeg,phpqdeg)/1.0d9	! dsigma/(dQ2 dxb dt dphi) in nb/GeV2--> ub/MeV2

	ntup%sigcm=sigma_eegamma*1.0d9
	
	jacobian_dis = (main%w**2-Mp**2)*xb*vertex%e%E/(2.0*pi*Mp*vertex%nu)
	
	peegamma=sigma_eegamma*jacobian*jacobian_dis

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

	real*8 xsec(46,56,36,21)

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

	   do i=1,36
	      phitab(i) = (i-1)*10.0
	   enddo

	   open(unit=10, file=filename, status='old')
	   read(10,'(a)') line !skip first line
	   
	   do i=1,46 ! Q2
	     do j=1,56		!xBj
		do k=1,36	!phipq
		   do l=1,21	!thpq
		      read(10,*) q2dum,xdum,tdum,phidum,thdum,xsec(i,j,k,l)
c		      write(6,*) q2dum,xdum,tdum,phidum,thdum,xsec(i,j,k,l)
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
		    do k=1,36	!phipq
		       if(phipq.ge.350.0) then
			  phimax=360.0
		       else
			  phimax=phitab(k+1)
		       endif
		       if(phipq.gt.phitab(k) .and. phipq.le.phimax) then
			  philo=phitab(k)
			  phihi=phitab(k+1)
			  dphi=phihi-philo
			  do l=1,21 !thpq
			     if(thpq.gt.thtab(l) .and. thpq.le.thtab(l+1)) then
				thlo=thtab(l)
				thhi=thtab(l+1)
				dth=thhi-thlo

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
			     endif ! theta
			  enddo ! theta
		       endif	! phi
		    enddo	! phi
		 endif		!xb
	      enddo		!xb
	   endif		! Q2
	enddo			! Q2

	dvcs_xsec=c

	return
	end
