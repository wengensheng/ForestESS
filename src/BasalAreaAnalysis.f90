 program growthana

 implicit none
 integer, parameter :: no_file=1
 integer, parameter :: items = 21
 integer, parameter :: planttraits = 18
 integer, parameter :: maxLayers =2
 integer, parameter :: maxPFTs  =3
 integer, parameter :: id_format = 103
 integer, parameter :: cohorts = 40
 integer, parameter :: maxyears=16009
 integer, parameter :: bins=21 ! 6
 integer, parameter :: bins2 = bins + 1
 real,    parameter :: PI = 3.1415926
 character(len=90) :: filein(no_file)
 character(len=60) :: filepath
 character(len=60) :: filein1,filein2,fileout1,fileout2
 character(len=60) :: comments
 integer,dimension(2,maxyears):: yearch
 real,   dimension(items,cohorts,maxyears):: dataarray,dataarray2
 real,   dimension(planttraits,maxPFTs,maxyears):: firstCC,L2firstCC
 real,   dimension(maxPFTs,maxyears):: BApft,dDBH,density
 real,   dimension(45,maxyears):: ecodata
 real,   dimension(maxyears):: HTstar,Layer1BM,Layer1Den,CAabv5m,totalBM,basalA
 real,   dimension(maxyears):: BMmort
 real,   dimension(planttraits,maxPFTs,maxLayers,maxyears):: meanvalues

 real,   dimension(bins2) :: DBHbins
 real,   dimension(bins,maxyears) :: DBHclasses,BMclasses

 integer :: totyears,m,n,i,j,k,iPFT,iLayer,layertype
 integer :: commentlines
 integer :: istat2,istat3,cc
 real :: DBH, GPP, NPP, plantC, soilC,plantN, soilN, mineralN, Nmin

 commentlines = 1 ! 91 ! 91
 filepath = 'output/' ! 'testruns/'
 filein1  = trim(filepath)//'Annual_cohorts.csv'
 filein2  = trim(filepath)//'EcosystemDynamics.csv'
 fileout1 = trim(filepath)//'Cohorts_mean.csv'

 open(211,file=filein1,status='old',ACTION='read',IOSTAT=istat2)

 open(13,file=fileout1)

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
       BMmort     = 0.
       basalA     = 0.
       BApft      = 0.0
       dDBH       = 0.
       density    = 0.
       m=0
       do
          m=m+1
          read(211,*,IOSTAT=istat3)yearch(1,m),yearch(2,m)
          !write(*,*)yearch(1,m),yearch(2,m)
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
              ! mean dDBH of the first layer
              if(dataarray(3,j,m)<2.0)then
                  dDBH(iPFT,m) = dDBH(iPFT,m) + dataarray(6,j,m) * dataarray(4,j,m)
                  density(iPFT,m) = density(iPFT,m) + dataarray(4,j,m)
              endif
              !if(dDBH(iPFT,m) == 0.0 .and. dataarray(3,j,m)<2.0)then
              !    dDBH(iPFT,m) = dataarray(6,j,m) * dataarray(4,j,m)
              !    density(iPFT,m) = density(iPFT,m) + dataarray(4,j,m)
              !endif
          enddo ! j, cohorts
          do i=1,maxPFTs
             if(density(i,m)>0.0) dDBH(i,m)=dDBH(i,m)/density(i,m)
          enddo
!         read the last cc and discard
          if(yearch(2,m) > cc)  &
              read(211,*,IOSTAT=istat3)(dataarray(i,cc+1,m),i=1,items)
          if(istat3<0)exit

          ! Output data
          !write(14,140)m,(BApft(i,m),i=1,3), &
          !               (dDBH(i,m), i=1,3), &
          !               (density(i,m),i=1,3)



       enddo ! End of calculating
       totyears=m-1
     else
           write(*,*)" open failed !"
     endif


! read in another file for combining data
   open(212,file=filein2,status='old',ACTION='read')
   read(212,*,IOSTAT=istat3)comments

! output to files
 fileout2=trim(filepath)//'WC-ResidentDC_Va1.5_SC33.csv'
 open(14,file=fileout2)
! write(14,'(360(a12,","))') 'year', &
!        'BA_DC','BA_EG','BA_fixer', &
!        'dD_DC','dD_EG','dD_fixer', &
!        'n_DC', 'n_EG','n_fixer'
 write(14,'(360(a12,","))') 'year', &
        'BA_DC','BA_EG','dD_DC','dD_EG', &
        'GPP', 'NPP','plantC','soilC',   &
        'plantN','soilN','mineralN','Nmin'

 do m = 1, totyears
!    write(13,103)m,HTstar(m),CAabv5m(m),Layer1Den(m),         &
!                   Layer1BM(m),BMmort(m),totalBM(m),          &
!                   basalA(m), (BApft(i,m),i=1,maxPFTs),        &
!                   (DBHclasses(i,m),i=1,bins),                &
!                   (BMclasses(i,m),i=1,bins),                 &
!                   (iPFT,(firstCC(i,iPFT,m),                  &
!                   i=1,planttraits),iPFT=1,maxPFTs),           &
!                   (iPFT,(L2firstCC(i,iPFT,m),                &
!                   i=1,planttraits),iPFT=1,maxPFTs)

     read(212,*,IOSTAT=istat3)n,(ecodata(i,m),i=1,32)
     GPP = ecodata(3,m)
     NPP = ecodata(4,m)
     plantC = ecodata(6,m)
     soilC  = ecodata(7,m)
     plantN = ecodata(8,m)
     soilN  = ecodata(9,m)
     mineralN = ecodata(28,m)
     Nmin   = ecodata(31,m)
     !if(m>500) &
     write(14,140)m,(BApft(i,m),i=1,2), (dDBH(i,m), i=1,2), &
                GPP, NPP, plantC, soilC,plantN, soilN, mineralN, Nmin
 enddo

140 format(1(I8,','),10(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
103 format(1(I8,','),11(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
102 format(1(I8,','),4(f15.4,','),6(f15.4,','), 2(I8,',',10(f15.4,',')), 12(2(I8,','),10(f15.4,',')) )
101 format(1(I8,','),4(f15.4,','),6(f15.4,','), 1(I8,',',10(f15.4,',')),  6(2(I8,','),10(f15.4,',')) )
close (12)
close (13)
close (212)
end
