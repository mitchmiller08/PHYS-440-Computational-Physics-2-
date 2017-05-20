	complex*16 a(600)
	integer i,k,inv,n
	open(unit=42,file="ddp.12.dat",status="unknown")
	open(unit=43,file="ddp.12.spec.dat",status="unknown")
	k=10
	n=2**k
	delT=0.4
	do i=1,n
	read(42,*) time,x
	a(i)=x
	end do
	inv=0
	call fft(a,k,inv)

	do j=1,n
	r1=real(a(j))
	ai=dimag(a(j))
	spec=sqrt(r1**2+ai**2)
	freq = j*2*3.141592653589793d0/60.000000000003780
	write(43,*) freq,spec
	end do
	stop
	end


	subroutine fft(a,m,inv)
	complex*16 a(600),u,w,t
	double precision ang,pi
	integer n,nd2,i,j,k,l,le,le1,ip
	parameter (pi=3.141592653589793d0)
	n=2**m
	nd2=n/2
	j=1
	do i=1,n-1
	if(i.lt.j) then
	t=a(j)
	a(j)=a(i)
	a(i)=t
	endif
	k=nd2
100	if(k.lt.j) then
	j=j-k
	k=k/2
	goto 100
	endif
	j=j+k
	end do
	le=1
	do l=1,m
	le1=le
	le=le+le
	u=(1.d0,0.d0)
	ang=pi/dble(le1)
	w=dcmplx(cos(ang),-sin(ang))
	if(inv.eq.1) w=dconjg(w)
	do j=1,le1
	do i=j,n,le
	ip=i+le1
	t=a(ip)*u
	a(ip)=a(i)-t
	a(i)=a(i)+t
	end do
	u=u*w
	end do
	end do
	if(inv.ne.1) then
	do i=1,n
	a(i)=a(i)/dble(n)
	end do
	endif
	end
