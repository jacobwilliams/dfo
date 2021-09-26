integer nvarmx, nconmx, nlinmx, nnlnmx, nfunmx
parameter (nvarmx=50, nconmx=50, nnlnmx=50, nlinmx=50, &
           nfunmx=10000)

!
!  CONSTRAINTS MODEL COEFFICIENTS
!
double precision ccon, lcon, qcon

common /mdlcon/  ccon(nconmx), lcon(nconmx*nvarmx), &
                 qcon(nconmx*nvarmx*nvarmx)
save /mdlcon/

double precision gmod, hmod, amat

common / mdlpar / gmod(nvarmx), hmod(nvarmx,nvarmx), &
                  amat(nlinmx, nvarmx)
save / mdlpar /

integer          useipopt, ncon, nnln, nlin 

common / mdldim / useipopt, ncon, nnln, nlin 
save / mdldim /
!
!  PARAMETERS NEED TO COMPUTE A MERIT FUNCTION AS OBJECTIVE FOR NPSOL
!

double precision conl, conu, penpar

logical          usemerit

common / merit / conl(nconmx), conu(nconmx), penpar, usemerit
save / merit /



