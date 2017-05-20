	dimension r(40000),e(40000),auto(40000)
	real lag_time, length
	double precision dseed,d2p31m,d2p31
	data d2p31m/2147483647.d0/
	data d2p31 /2147483711.d0/
	open(unit=42,file='range1.2.dat',status='unknown')
	open(unit=43,file='range2.2.dat',status='unknown')
	open(unit=44,file='range3.2.dat',status='unknown')
	
	delta_time=0.00325
	length=6032
	alpha=0.05
	time=0.0
	w0=1.5
	c=0.05

	do i=1,35201
	time=time+delta_time
	t=time
	if(i.ge.length) then
	r(i)=0.0
	else
	r(i)=sin((w0+c*t)*t)
	endif
	write(42,*) time,r(i)
	end do

	time=0.0
	do i=1,35201
	time=time+delta_time
	data dseed/3156789.d0/
	dseed=mod(16807.d0*dseed,d2p31m)
	rannos=dseed/d2p31
	if(i.le.14000) e(i)=rannos-0.5
	if(i.ge.14000+length) e(i)=rannos-0.5
	if(i.gt.14001.and.i.lt.14000+length) e(i)=alpha*r(i-14001)+rannos-0.5
	if(i.gt.24001.and.i.lt.24000+length) e(i)=alpha*r(i-24001)+rannos-0.5
	write(43,*) time, e(i)
	end do

	do j=1,32000
	lag_time=j*delta_time
	
	sum=0.0
	do i=1,length
	sum=sum+r(i)*e(i+j)*delta_time
	end do
	duration=length*delta_time
	auto(j)=sum/duration

	write(44,*) lag_time,auto(j)
	end do

	stop
	end
