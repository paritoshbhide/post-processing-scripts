program distribution

!module Particles
implicit none

type Particle_Data
 real*8::Cx
 real*8::Cy
 real*8::Cz
 real*8::x_P
 real*8::y_P
 real*8::z_P
 integer::active
end type Particle_Data

type(Particle_Data),dimension(500000)::Particle
!end module Particles


!******************************************
!All constants here
!******************************************
module constants
!implicit none
real*8, parameter :: avogadro = 6.02214129d+23
real*8, parameter :: Pi=3.1415926535897
real*8, parameter :: Kb=1.3806488d-23
real*8, parameter :: hplanck=6.62606957d-34
end module constants

REAL:: random_numb

integer:: i,j,k,particle_index,N_Particle,delete
real*8 :: Temperature,mass_dia,R_gas_constant

real*8 :: u_bulk,v_bulk,w_bulk, C_peculiar1, C_Peculiar2

real*8,dimension(1000)::Maxwellfx,Maxwellfy,Maxwellfz,Maxwell_C
real*8::A1,A2
en

!*******************************************************
!Binning Width And Parameter
type bin
integer:: np_bin
real*8 :: bin_v
end type bin

type(bin), allocatable :: bin_vx, bin_vy, bin_vz
!**********************************************************





character(len=25)::Filename_in
character(len=5)::fmt
character(8)::numbr

write(*,*) 'Enter Number of Particles'
read(*,*), N_particle

write(*,*) 'Enter file name'
read(*,*), Filename_in

R_gas_Constant=208
bin_velocity%bin_Number=100



open(unit=2, file=Filename_in)
do particle_index=1,N_particle
     read(2,*) Particle(particle_index)%Cx,Particle(particle_index)%Cy,Particle(particle_index)%Cz
enddo
close(2)


!Taking the average
u_bulk=0.0d0
v_bulk=0.0d0
w_bulk=0.0d0
do particle_index=1,N_Particle
   u_bulk=u_bulk+Particle(particle_index)%Cx
   v_bulk=v_bulk+Particle(particle_index)%Cy
   w_bulk=w_bulk+Particle(particle_index)%Cz
end do

u_bulk=u_bulk/N_Particle
v_bulk=v_bulk/N_Particle
w_bulk=w_bulk/N_Particle


print*, u_bulk,' ',v_bulk,' ',w_bulk



!Calculating temperature
Temperature=0.0d0
do particle_index=1,N_particle
    Temperature=(Particle(particle_index)%Cx-u_bulk)**2 +&
                (Particle(particle_index)%Cy-v_bulk)**2 +&
                (Particle(particle_index)%Cz-w_bulk)**2 + Temperature
end do

Temperature=Temperature/(3.0d0*R_gas_Constant*N_Particle)

print*,temperature
!Binning of All Velocities
bin_velocity%bin_width=(maxval(Particle(:)%Cx)-minval(Particle(:)%Cx))/100
bin_velocity%bins_C(1)=minval(Particle(:)%Cx)

do i=1,100
bin_velocity%bins_C(i+1)=bin_velocity%bins_C(i)+bin_velocity%bin_width
end do


bin_velocity%bins_Cx(:)=0
bin_velocity%bins_Cy(:)=0
bin_velocity%bins_Cy(:)=0

do i=1,99
       do particle_index=1,N_Particle
           if(Particle(particle_index)%Cx .ge. bin_velocity%bins_C(i)) then
           if(Particle(particle_index)%Cx .lt. bin_velocity%bins_C(i+1)) then
                bin_velocity%bins_Cx(i)=bin_velocity%bins_Cx(i)+1
         end if
         end if
    end do
end do

!Generating Maxwellian From Analytical Equation
u_bulk=1000.0d0
v_bulk=0.0d0
w_bulk=0.0d0
mass_dia=39.0d0*1d-3/avagadro
Temperature=1000.0d0
A1=sqrt(mass_dia/(2.0d0*Pi*kb*Temperature))
do i=1,99
   !Maxwell_C(i)=bin_velocity%bins_C(1)+bin_velocity%bin_width/10.0d0
   Maxwellfx(i)=A1*Exp(-(mass_dia/(2.0d0*kb*temperature))*(bin_velocity%bins_C(i)-u_bulk)**2.0d0)
   Maxwellfy(i)=A1*Exp(-(mass_dia/(2.0d0*kb*temperature))*(bin_velocity%bins_C(i)-v_bulk)**2.0d0)
   Maxwellfz(i)=A1*Exp(-(mass_dia/(2.0d0*kb*temperature))*(bin_velocity%bins_C(i)-w_bulk)**2.0d0)
end do

delete=0
dummy=0
do i=1,100
particle_distro(i)=bin_velocity%bins_Cx(i)
dummy=(bin_velocity%bin_width)*bin_velocity%bins_Cx(i)
delete=delete+dummy
end do



print*,delete

particle_distro(:)=particle_distro(:)/delete
!bin_velocity%bins_Cx(:)=bin_velocity%bins_Cx(:)/delete

fmt='(I5)'
write(numbr,fmt) 10
filename="distroT1000u1000P1.dat"
open(unit=1, file=filename, status='replace',action="write")
3 format(f20.15,a,f20.15)
write(1,*)"variables=C,Particles,Maxwell"
do i=1,100
    write(1,*),bin_velocity%bins_C(i),particle_distro(i),Maxwellfx(i)
end do
close(1)

filename="Maxwell.dat"
open(unit=2, file=filename, status='replace',action="write")
2 format(f20.15,a,f20.15)
write(2,*)"variables=C,Cx"
do i=1,100
    write(2,*),bin_velocity%bins_C(i),Maxwellfx(i)
end do
close(2)

end program distribution
