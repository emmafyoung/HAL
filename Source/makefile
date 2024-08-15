SHELL=/bin/csh

# PROGRAM NAME
main=ibm
#MODULE NAME
mod1=param
#EXTERNAL SUBROUTINE NAMES
sr1=Filope
sr4=Diapause_depth
sr5=Check
sr6=Mod_param
sr7=Input_bio
sr8=Ran
sr9=Setup_output_netcdf
sr10=Input_phys
sr11=Input_NEMO_netcdf
sr12=Input_POM
sr13=Scoordset
sr14=Sigmaset
sr15=Output_netcdf
sr16=Input_bio2
sr17=Convert
sr18=Read_POLCOMS
sr19=Read_OCCAM
sr20=Read_POM
sr21=Read_NEMO_netcdf
sr22=Day_light
sr23=Bio_Calc_Egg
sr24=Bio_Calc_Larvae
sr25=Bio_Calc_Krill
sr26=Bio_Calc_Copepod
sr27=Calc_dens

FC=ifort
CO=-O2

objects1 = $(mod1).o $(sr5).o $(sr8).o $(sr1).o $(sr4).o $(sr6).o $(sr7).o $(sr9).o $(sr10).o $(sr11).o $(sr12).o $(sr13).o $(sr14).o $(sr15).o $(sr16).o $(sr17).o $(sr18).o $(sr19).o $(sr20).o $(sr21).o $(sr22).o $(sr23).o $(sr24).o $(sr25).o $(sr26).o $(sr27).o $(main).o

$(main):	$(objects1)
	$(FC) -o $(main) $(CO) $(objects1) -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(mod1).o:	$(mod1).F
	$(FC) -c $(CO) $(mod1).F
$(sr5).o:	$(sr5).F
	$(FC) -c $(CO) $(sr5).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(sr8).o:	$(sr8).F
	$(FC) -c $(CO) $(sr8).F
$(sr1).o:	$(sr1).F
	$(FC) -c $(CO) $(sr1).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(sr4).o:	$(sr4).F
	$(FC) -c $(CO) $(sr4).F -shared-intel
$(sr6).o:	$(sr6).F
	$(FC) -c $(CO) $(sr6).F
$(sr7).o:	$(sr7).F
	$(FC) -c $(CO) $(sr7).F
$(sr9).o:	$(sr9).F
	$(FC) -c $(CO) $(sr9).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(sr10).o:	$(sr10).F
	$(FC) -c $(CO) $(sr10).F
$(sr11).o:	$(sr11).F
	$(FC) -c $(CO) $(sr11).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(sr12).o:	$(sr12).F
	$(FC) -c $(CO) $(sr12).F
$(sr13).o:	$(sr13).F
	$(FC) -c $(CO) $(sr13).F
$(sr14).o:	$(sr14).F
	$(FC) -c $(CO) $(sr14).F
$(sr15).o:	$(sr15).F
	$(FC) -c $(CO) $(sr15).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(sr16).o:	$(sr16).F
	$(FC) -c $(CO) $(sr16).F
$(sr17).o:	$(sr17).F
	$(FC) -c $(CO) $(sr17).F
$(sr18).o:	$(sr18).F
	$(FC) -c $(CO) $(sr18).F
$(sr19).o:	$(sr19).F
	$(FC) -c $(CO) $(sr19).F
$(sr20).o:	$(sr20).F
	$(FC) -c $(CO) $(sr20).F
$(sr21).o:	$(sr21).F
	$(FC) -c $(CO) $(sr21).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel
$(sr22).o:	$(sr22).F
	$(FC) -c $(CO) $(sr22).F
$(sr23).o:	$(sr23).F
	$(FC) -c $(CO) $(sr23).F
$(sr24).o:	$(sr24).F
	$(FC) -c $(CO) $(sr24).F
$(sr25).o:	$(sr25).F
	$(FC) -c $(CO) $(sr25).F
$(sr26).o:	$(sr26).F
	$(FC) -c $(CO) $(sr26).F
$(sr27).o:	$(sr27).F
	$(FC) -c $(CO) $(sr27).F
$(main).o:	$(main).F
	$(FC) -c $(CO) $(main).F -L$(NETCDFDIR) -lnetcdf -lnetcdff -I$(NETCDFINCLUDE) -shared-intel

