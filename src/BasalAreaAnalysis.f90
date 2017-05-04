 program growthana

 implicit none
 integer, parameter :: no_file=1
 integer, parameter :: items = 22
 integer, parameter :: planttraits = 19
 integer, parameter :: maxLayers =2
 integer, parameter :: maxPFTs  =3
 integer, parameter :: id_format = 103
 integer, parameter :: cohorts = 90
 integer, parameter :: maxyears=8001
 integer, parameter :: bins=21 ! 6
 integer, parameter :: bins2 = bins + 1 
 real,    parameter :: PI = 3.1415926
 character(len=90) :: filein(no_file)
 character(len=60) :: field(no_file),fileout1,fileout2,filepath
 character(len=60) :: comments
 integer,dimension(2,maxyears):: yearch
 real,   dimension(items,cohorts,maxyears):: dataarray,dataarray2
 real,   dimension(planttraits,maxPFTs,maxyears):: firstCC,L2firstCC
 real,   dimension(maxPFTs,maxyears):: BApft
 real,   dimension(maxyears):: HTstar,Layer1BM,Layer1Den,CAabv5m,totalBM,basalA
 real,   dimension(maxyears):: BMmort
 real,   dimension(planttraits,maxPFTs,maxLayers,maxyears):: meanvalues 

 real,   dimension(bins2) :: DBHbins
 real,   dimension(bins,maxyears) :: DBHclasses,BMclasses

 integer :: totyears,m,n,i,j,k,iPFT,iLayer,layertype
 integer :: commentlines
 integer :: istat2,istat3,cc
 real :: DBH

 commentlines = 1 ! 91 ! 91
 filepath = '' ! 'testruns/'
 field(1) = 'Annual_cohorts' !

! fileout1=trim(field(1))//'_mean.csv'
 fileout2=trim(field(1))//'_basalarea.csv'

 open(14,file=fileout2)
 write(14,'(360(a12,","))') 'year','Batot','Deciduous','Evergreen'

! file names
 do i=1,no_file
     filein(i)=trim(filepath)//trim(field(i))//'.txt'  !trim(filepath)//trim(field(i))
 enddo

! initiate bins
DBHbins(1) = 0.0 ! 0.05
do i=2,bins
    DBHbins(i) = DBHbins(i-1)+0.1 !0.1
enddo
DBHbins(1) = 0.025
DBHbins(bins+1) = 9999.
DBHclasses = 0.0
BMclasses = 0.0   
! read in file
   open(211,file=filein(1),status='old',ACTION='read',IOSTAT=istat2)
   if(istat2==0)then
       write(*,*)" open successfully!"
       do i=1,commentlines !91
          read(211,*,IOSTAT=istat3)comments
       enddo

       firstCC = -9999.
       L2firstCC = -9999.
       meanvalues = 0.
       HTstar     = 0. 
       Layer1BM   = 0.
       Layer1Den  = 0.
       CAabv5m    = 0.
       totalBM    = 0.
       BMmort      = 0.
       basalA     = 0.
       BApft      = 0.0
       m=0
       do 
          m=m+1
          read(211,*,IOSTAT=istat3)yearch(1,m),yearch(2,m)
          write(*,*)yearch(1,m),yearch(2,m)
          if(istat3<0)exit
          cc = yearch(2,m)

          if(yearch(2,m) > 4) cc = yearch(2,m) -1 ! skip the last cohort

          do j=1,cc
              read(211,*,IOSTAT=istat3)(dataarray(i,j,m),i=1,items)
              if(istat3<0)exit
!              dataarray(4,j,m)=dataarray(4,j,m)/10000.

!             CA of plants above 5 meters
              if(dataarray(8,j,m) > 0.5) then ! height
                  CAabv5m(m) = CAabv5m(m) + dataarray(9,j,m) ! crown area
              endif
!             first layer boimass and HT*
              if(dataarray(3,j,m) < 2.) then
                  HTstar(m)    = dataarray(7,j,m)
                  Layer1Den(m) = Layer1Den(m) + dataarray(4,j,m)
                  Layer1BM(m)  = Layer1BM(m)  + dataarray(4,j,m) * &
                                (dataarray(10,j,m)+dataarray(11,j,m))/10000
              endif
!             total biomass
              totalBM(m)  = totalBM(m)  + (dataarray(11,j,m)+dataarray(10,j,m))  &
                            * dataarray(4,j,m)/10000
!             woody biomass residence time
              DBH = dataarray(7,j,m)
              iLayer=dataarray(3,j,m)
              iPFT  =dataarray(2,j,m)-1


!             calculate means or sums
              dataarray2(:,j,m) = dataarray(:,j,m)
              dataarray2(6:items,j,m) = dataarray(6:items,j,m) * dataarray(4,j,m)

!             size classes distribution, including total individuals; biomass distribution
              do i=1,bins
                  if(dataarray(7,j,m)>DBHbins(i) .and. dataarray(7,j,m)<DBHbins(i+1) ) then
                     DBHclasses(i,m) = DBHclasses(i,m) + dataarray(4,j,m)
                     BMclasses(i,m)  = BMclasses(i,m)  + dataarray(4,j,m)*         &
                                       (dataarray(11,j,m)+dataarray(10,j,m)) /10000
                   endif
              enddo

              iLayer=dataarray(3,j,m)
              iPFT  =dataarray(2,j,m)-1  ! dataarray(2,j,m)-2
!             for all the variables other than dDBH, summation ! from nindivs
!              meanvalues(:,iPFT,iLayer,m)=meanvalues(:,iPFT,iLayer,m)+dataarray2(4:items,j,m)


!             Basal area with DBH > 0.1 m
              if(dataarray(7,j,m)> 0.0)then
                  basalA(m)   = basalA(m)   + 0.25*PI*dataarray(7,j,m)**2 * dataarray(4,j,m)
                  BApft(iPFT,m) = BApft(iPFT,m) + 0.25*PI*dataarray(7,j,m)**2 * dataarray(4,j,m)
              endif

          enddo
!         read the last cc and discard
          if(yearch(2,m) > cc)  &
              read(211,*,IOSTAT=istat3)(dataarray(i,cc+1,m),i=1,items)
          if(istat3<0)exit

          ! Output data
          write(14,140)m,basalA(m),BApft(2,m),BApft(3,m)

        enddo ! End of calculating
        totyears=m-1
     else
           write(*,*)" open failed !"
     endif

140 format(1(I8,','),10(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
103 format(1(I8,','),11(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
102 format(1(I8,','),4(f15.4,','),6(f15.4,','), 2(I8,',',10(f15.4,',')), 12(2(I8,','),10(f15.4,',')) )
101 format(1(I8,','),4(f15.4,','),6(f15.4,','), 1(I8,',',10(f15.4,',')),  6(2(I8,','),10(f15.4,',')) )

close (211)
close (14)
end
