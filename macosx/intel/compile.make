#
#
#  ###################################################################
#  ##                                                               ##
#  ##  compile.make  --  compile all the TINKER modules for OpenMP  ##
#  ##         (Intel Fortran Compiler for Mac OSX Version)          ##
#  ##                                                               ##
#  ###################################################################
#
#
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp active.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp alchemy.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp analysis.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp analyze.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp angles.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp anneal.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp archive.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp attach.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp bar.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp basefile.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp beeman.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp bicubic.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp bitors.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp bonds.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp born.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp bounds.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp bussi.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp calendar.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp center.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp chkpole.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp chkring.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp chkxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp cholesky.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp clock.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp cluster.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp column.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp command.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp connect.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp connolly.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp control.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp correlate.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp crystal.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp cspline.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp cutoffs.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp deflate.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp delete.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp diagq.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp diffeq.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp diffuse.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp distgeom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp document.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp dynamic.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangang.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangang1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangang2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangang3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangle.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangle1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangle2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eangle3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebond.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebond1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebond2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebond3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebuck.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebuck1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebuck2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ebuck3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echarge.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echarge1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echarge2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echarge3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echgdpl.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echgdpl1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echgdpl2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp echgdpl3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp edipole.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp edipole1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp edipole2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp edipole3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egauss.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egauss1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egauss2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egauss3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egeom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egeom1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egeom2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp egeom3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ehal.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ehal1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ehal2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ehal3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimprop.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimprop1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimprop2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimprop3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimptor.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimptor1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimptor2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eimptor3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elj.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elj1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elj2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elj3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp embed.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emetal.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emetal1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emetal2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emetal3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emm3hb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emm3hb1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emm3hb2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp emm3hb3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp empole.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp empole1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp empole2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp empole3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp energy.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopbend.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopbend1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopbend2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopbend3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopdist.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopdist1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopdist2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eopdist3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp epitors.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp epitors1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp epitors2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp epitors3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp erf.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp erxnfld.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp erxnfld1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp erxnfld2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp erxnfld3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esolv.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esolv1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esolv2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esolv3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrbnd.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrbnd1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrbnd2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrbnd3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrtor.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrtor1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrtor2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp estrtor3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etors.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etors1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etors2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etors3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etortor.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etortor1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etortor2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp etortor3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eurey.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eurey1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eurey2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp eurey3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp evcorr.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp extra.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp extra1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp extra2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp extra3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp fatal.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp fft3d.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp fftpack.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp field.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp final.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp flatten.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp freeunit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp gda.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp geometry.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getint.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getkey.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getmol.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getmol2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getnumb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getpdb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getprm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getref.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getstring.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp gettext.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getword.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp getxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ghmcstep.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp gradient.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp gradrgd.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp gradrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp groups.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp grpline.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp gyrate.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp hessian.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp hessrgd.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp hessrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp hybrid.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp image.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp impose.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp induce.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp inertia.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp initatom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp initial.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp initprm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp initres.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp initrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp insert.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp intedit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp intxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp invbeta.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp invert.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp jacobi.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kangang.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kangle.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp katom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kbond.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kcharge.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kdipole.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kewald.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kextra.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kgeom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kimprop.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kimptor.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kinetic.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kmetal.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kmpole.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kopbend.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kopdist.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp korbit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kpitors.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kpolar.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ksolv.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kstrbnd.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kstrtor.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ktors.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ktortor.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kurey.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp kvdw.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp lattice.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp lbfgs.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp lights.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp makeint.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp makeref.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp makexyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp maxwell.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp mdinit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp mdrest.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp mdsave.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp mdstat.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp mechanic.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp merge.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp minimize.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp minirot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp minrigid.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp molecule.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp molxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp moments.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp monte.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp mutate.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp nblist.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp newton.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp newtrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp nextarg.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp nexttext.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp nose.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp nspline.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp nucleic.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp number.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp numeral.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp numgrad.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp ocvm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp openend.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp optimize.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp optirot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp optrigid.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp optsave.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp orbital.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp orient.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp orthog.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp overlap.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp path.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pdbxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp picalc.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pmestuff.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pmpb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp polarize.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp poledit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp polymer.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp potential.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp precise.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pressure.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prmedit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prmkey.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp promo.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp protein.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtdyn.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prterr.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtint.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtmol2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtpdb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtprm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtseq.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp prtxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pss.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pssrigid.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp pssrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp quatfit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp radial.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp random.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp rattle.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readdyn.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readgau.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readint.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readmol.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readmol2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readpdb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readprm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readseq.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp readxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp replica.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp respa.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp rgdstep.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp rings.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp rmsfit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp rotlist.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp rotpole.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp saddle.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp scan.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp sdstep.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp search.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp server.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp shakeup.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp sigmoid.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp sktstuff.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp sniffer.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp sort.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp spacefill.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp spectrum.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp square.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp suffix.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp superpose.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp surface.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp surfatom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp switch.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp sybylxyz.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp temper.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp testgrad.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp testhess.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp testpair.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp testpol.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp testrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp timer.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp timerot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp tncg.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp torphase.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp torque.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp torsfit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp torsions.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp trimtext.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp unitcell.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp valence.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp verlet.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp version.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp vibbig.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp vibrate.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp vibrot.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp volume.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xtalfit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xtalmin.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xyzatm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xyzedit.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xyzint.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xyzpdb.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp xyzsybyl.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp zatom.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp klfmm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elfmm.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elfmm1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elfmm2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp elfmm3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esnbnd.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esnbnd1.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esnbnd2.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp esnbnd3.f
ifort -c -O3 -axSSSE3 -no-ipo -no-prec-div -mdynamic-no-pic -w -assume cc_omp -openmp iccglfse.f
