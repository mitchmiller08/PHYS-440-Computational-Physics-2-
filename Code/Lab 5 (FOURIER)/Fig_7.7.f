	Double Precision psi3(2048), deltat
	Complex*16 ExpT(2048), ExpV(2048), im, phi(2048), psi(2048)
	Complex*16 psi2(2048)
	Double precision V(2048), KE, ko, dt, maxtime, k, Energy
	Double precision hbar, hbar2, deltak, time, pi, x(2048), h
	integer q
	
	open(unit=21,file='wave_new.7.dat',status='unknown')
	open(unit=22,file='wave_old.7.dat',status='unknown')
	open(unit=23,file='potent_new.7.dat',status='unknown')

	im=(0,1)
	pi=3.141592653589793d0
	hbar=6.5821220d0
	hbar2=7.6199682d0
	mass=1.d0
	Energy=0.d0
	ko=sqrt((2*mass*Energy)/hbar2)
	deltat=5.0d-2
	dt=5.0d-2
c	Collision around 37*dt
	maxtime=45*dt
	p=7
	N=128
	h=0.3125d0
	deltak=2*pi/(N*h)
	
	x(1)=-20
	DO q=1,N
c	May need to change 4.0d0 to 4.0d0*h*h
	psi(q)=0.6316187777*cdexp((im*ko*(x(q)))-((x(q))**2)/(4.0d0))
	x(q+1)=x(q)+h
	wavemag=dsqrt((dreal(psi(q)))**2+(dimag(psi(q)))**2)
	write(22,*) x(q),wavemag
	END DO

	DO q=1,N/2
		k=dble(q-1)*deltak
		KE=k*k*hbar2/(2.d0*mass)
		ExpT(q)=cdexp(-im*KE*deltat/hbar)
		k=-dble(q)*deltak
		KE=k*k*hbar2/(2.d0*mass)
		ExpT(N+1-q)=cdexp(-im*KE*deltat/hbar)
	END DO

	Do q=1,128
		V(q)=0.0
	END Do

	do q=1,N
		write(23,*) x(q), V(q)/200.d0
		ExpV(q)=cdexp(-im*V(q)*deltat/(2*hbar))
c		write(23,*) x(q), ExpV(q)
	end do

	time=0.d0
	inc=0
1000	time=time+dt
	inc=inc+1

	do q=1,N
		phi(q)=ExpV(q)*psi(q)
	end do

	inv=0
	call FFT(phi,p,inv)

	do q=1,N
		phi(q)=ExpT(q)*phi(q)
	end do

	inv=1
	call fft(phi,p,inv)

	do q=1,N
		psi(q)=ExpV(q)*phi(q)
	end do

c	if(inc.EQ.11) then
c		do q=1,N
c			psi2(q)=Dconjg(psi(q))
c			psi3(q)=psi(q)*psi2(q)
c			write(21,*) x(q), psi3(q)
c		end do
c	end if

	if(time.lt.maxtime) goto 1000

c	if(inc.EQ.1000) then
		do q=1,N
			psi2(q)=Dconjg(psi(q))
			psi3(q)=psi(q)*psi2(q)
			write(21,*) x(q), psi3(q)
		end do
c	end if

	stop
	end

	Subroutine FFT(a,m,inv)
	Complex*16 a(5000),u,w,t
	double precision ang, pi
	Integer N,Nd2,i,j,r,l,le,le1,ip
	Parameter (pi=3.141592653589793d0)

	N=128
	Nd2=N/2
	j=1
	do i=1,N-1
		if(i.lt.j) then
			t=a(j)
			a(j)=a(i)
			a(i)=t
		end if
		r=Nd2
100		if(r.lt.j) then
			j=j-r
			r=r/2
			goto 100
		end if
		j=r+j
	end do

	le=1
	do l=1,7
		le1=le
		le=le+le
		u=(1.d0,0.d0)
		ang=pi/dble(le1)
		w=Dcmplx(cos(ang),-sin(ang))
		if(inv.eq.1) w=Dconjg(w)

		do j=1,le1
			do i=j,N,le
				ip=i+le1
				t=a(ip)*u
				a(ip)=a(i)-t
				a(i)=a(i)+t
			end do
			u=u*w
		end do

	end do

	if(inv.ne.1) then
		do i=1,N
			a(i)=a(i)/dble(N)
		end do
	end if
	end
