      INTEGER NVARMX, NCONMX, NLINMX, NNLNMX, NFUNMX
      PARAMETER (NVARMX=50, NCONMX=50, NNLNMX=50, NLINMX=50, 
     +           NFUNMX=10000)

C
C  CONSTRAINTS MODEL COEFFICIENTS
C
      DOUBLE PRECISION CCON, LCON, QCON
    
      COMMON /MDLCON/  CCON(NCONMX), LCON(NCONMX*NVARMX),
     .                 QCON(NCONMX*NVARMX*NVARMX)
      SAVE /MDLCON/

      DOUBLE PRECISION GMOD, HMOD, AMAT
 
      COMMON / MDLPAR / GMOD(NVARMX), HMOD(NVARMX,NVARMX), 
     .                  AMAT(NLINMX, NVARMX)
      SAVE / MDLPAR /
  
      INTEGER          USEIPOPT, NCON, NNLN, NLIN 
      
      COMMON / MDLDIM / USEIPOPT, NCON, NNLN, NLIN 
      SAVE / MDLDIM /
C
C  PARAMETERS NEED TO COMPUTE A MERIT FUNCTION AS OBJECTIVE FOR NPSOL
C

      DOUBLE PRECISION CONL, CONU, PENPAR
 
      LOGICAL          USEMERIT

      COMMON / MERIT / CONL(NCONMX), CONU(NCONMX), PENPAR, USEMERIT
      SAVE / MERIT /


